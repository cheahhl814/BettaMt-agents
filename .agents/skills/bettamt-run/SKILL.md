---
name: bettamt-run
description: Run BettaMt end-to-end. Use when a user wants to assemble a mitogenome and may not know the full workflow. Walks them through preflight → nextflow run → qc/debug. Triggers: "run BettaMt", "assemble my mitogenome", "I have reads and a reference", "bettamt run", "start a BettaMt run", "bettamt pipeline".
---

# bettamt-run — BettaMt end-to-end orchestrator

This skill is a thin router. It does not duplicate logic from the other three `bettamt-*` skills. Its job is to ask: **"what stage is the user at, and which skill should they invoke next?"**

## 0. Resolve the pipeline and run directory

Same convention as the other skills: `$BETA_MT_HOME` env var, or sibling `../BettaMt`.

```bash
: "${BETA_MT_HOME:=$(cd "$(dirname "$SKILL_DIR")/../BettaMt" 2>/dev/null && pwd)}"
test -f "$BETA_MT_HOME/main.nf" || { echo "Cannot locate BettaMt. Set BETA_MT_HOME."; exit 1; }

# Run dir is wherever the user wants outputs to land. Default to the pipeline root.
: "${RUN_DIR:=$BETA_MT_HOME}"
```

## 1. Determine the user's stage

Try to detect automatically **before** asking:

```bash
# A polished mitogenome exists → run succeeded, time for QC
if ls "$RUN_DIR"/results/polish/*_racon.fasta  "$RUN_DIR"/results/polish/*_polca.fasta >/dev/null 2>&1; then
    STAGE="qc"
# Annotation output exists and has been checked → done
elif [ -f "$RUN_DIR/annotation-report.md" ]; then
    STAGE="done"
# Annotation output exists but no annotation report yet → annotate-qc
elif [ -f "$RUN_DIR"/results/annotation/annotated_genes.fasta ]; then
    STAGE="annotate-qc"
# A diagnosis already exists → debug done, walk through it
elif [ -f "$RUN_DIR/diagnosis.md" ]; then
    STAGE="review-diagnosis"
# Params ready → ready to run
elif [ -f "$RUN_DIR/params.json" ]; then
    STAGE="run"
# Nothing yet → preflight
else
    STAGE="preflight"
fi
```

If auto-detection is ambiguous, ask one short question:

> Are you starting a new run, or continuing a previous one?
> 
> - new run (no params.json)
> - I have params.json already
> - I just ran it and it failed
> - I just ran it successfully

## 2. Route to the right skill

| Stage              | Action                                                                                                   |
| ------------------ | -------------------------------------------------------------------------------------------------------- |
| `preflight`        | Invoke `bettamt-preflight`. It produces `params.json` and the run command. Come back here when done.     |
| `run`              | Show the run command (see §3). Get explicit confirmation. Then run.                                      |
| `debug`            | Invoke `bettamt-debug`. It writes `diagnosis.md`. Come back when done.                                   |
| `qc`               | Invoke `bettamt-qc`. It writes `report.md`.                                                              |
| `annotate-qc`      | Invoke `bettamt-annotate-qc`. It writes `annotation-report.md`.                                          |
| `review-diagnosis` | Read `$RUN_DIR/diagnosis.md` with the user, recommend Option A from its "Recommended next action" block. |

**Do not skip preflight.** Even if the user supplies params manually, recommend running `bettamt-preflight` first to validate the choices and produce the audit trail.

## 3. The run command (only fires at the `run` stage)

### Show the user the command

```bash
cd "$BETA_MT_HOME"

# Local (workstation):
nextflow run main.nf -profile local  -params-file "$RUN_DIR/params.json"

# SLURM (cluster):
nextflow run main.nf -profile slurm  -params-file "$RUN_DIR/params.json"
```

Ask: **"This will start a Nextflow run with the params above. OK to proceed? (yes / no / dry-run)"**

### Recommend dry-run for first-time users

```bash
nextflow run main.nf -profile local -params-file "$RUN_DIR/params.json" -dry-run
```

The `-dry-run` flag shows what Nextflow *would* execute (the DAG, with all process directives resolved) without actually running anything. Invaluable for inexperienced users to see the planned graph, catch typos in `params.json`, and understand resource allocation before spending SLURM credits.

### Resume is a separate, distinct command

If the user says "I already ran this, just resume", do NOT use the bare `nextflow run` command. Use:

```bash
nextflow run main.nf -profile slurm -params-file "$RUN_DIR/params.json" -resume
```

`-resume` picks up cached steps from `.nextflow/` and re-runs only what failed or what changed. Always suggest it as the *first* retry after a failure.

## 4. After the run

| Outcome                          | Next skill                                                        |
| -------------------------------- | ----------------------------------------------------------------- |
| Exit 0, polished FASTA exists    | `bettamt-qc`                                                      |
| Exit 0, annotation files exist   | `bettamt-annotate-qc` (if `--ref_gff`/`--ref_gb` was supplied)    |
| Exit non-zero, no `diagnosis.md` | `bettamt-debug`                                                   |
| `diagnosis.md` already exists    | walk the user through Option A from it                            |
| User says "looks good" after QC  | done — point them to `$RUN_DIR/report.md` and `$RUN_DIR/results/` |

If the run succeeded and produced both polished FASTA and annotated genes, the recommended order is:
**`bettamt-qc` first** (checks assembly integrity), then **`bettamt-annotate-qc`** (checks gene transfer). Annotation QC is meaningless on a broken assembly, so do assembly QC first.

After invoking any next-step skill, **explicitly hand off** — don't make the user discover the next action. Example:

> I've handed off to `bettamt-qc`. It will read the polished FASTA and the read evidence, then write `report.md` to your run directory. Say the word and I'll invoke it.

## 5. Important rules

- **Never** run the pipeline without explicit user confirmation
- **Never** edit `params.json` after preflight — if values are wrong, re-run preflight
- **Always** recommend `-dry-run` for first-time users (unless they explicitly say "just run it")
- **Always** offer to invoke the next skill — don't leave the user to find it
- **Never** run two nextflow runs in the same `$RUN_DIR` concurrently — `-resume` cannot disambiguate them. If the user wants to start over, move `.nextflow/` aside first or use a new run directory.
- **If `$BETA_MT_HOME` doesn't resolve, ask** before guessing. A wrong pipeline home = wasted compute.

## 6. Output contract — no artifact

This skill produces **no file** of its own. Its only outputs are:

- The run command (printed to the user)
- The hand-off message naming the next skill

All artifacts (`params.json`, `diagnosis.md`, `report.md`) are produced by the other three skills. Don't try to write them from here.

## 7. Common follow-ups

| User says                    | What to do                                                                                                                               |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| "Looks good" after QC        | Summarize: pipeline version, polished FASTA path, key QC verdicts                                                                        |
| "Looks good" after annotate-qc | Summarize: gene count, gene types, any BLAST warnings; point to `annotation-report.md`                                                  |
| "Can I publish this?"        | Point them to `$RUN_DIR/report.md`, `$RUN_DIR/annotation-report.md`, and `$RUN_DIR/params.rationale.md` — all are audit-trail artifacts suitable for supplementary material |
| "I want to run more samples" | Recommend making a new `$RUN_DIR` per sample, symlink `$BETA_MT_HOME`                                                                    |
| "How do I cite this?"        | `cat "$BETA_MT_HOME/README.md"` — citations block at the bottom                                                                          |
