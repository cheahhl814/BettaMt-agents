---
name: bettamt-run
description: Phase 2 of the bettamt skill bundle — execute the BettaMt mitogenome workflow using the runner selected at preflight. Default is the manual bash-recipes path (runners/bash-recipes/run_ont.sh / run_illumina.sh, per-stage headers, pinned pixi env); the opt-in Nextflow runner (runners/nextflow-runner/main.nf with -profile/-c override) is used for cluster or resume-heavy workflows. Also routes between pipeline stages, tracks completion state from the filesystem, and orchestrates annotate_genes.sh when a reference annotation was supplied. Triggers "run BettaMt", "bettamt run", "execute mitogenome assembly", "assemble this mito genome".
version: 1.0.0
updated: "2026-10-18"
triggers:
  - "run BettaMt"
  - "bettamt run"
  - "execute mitogenome assembly"
  - "assemble this mito genome"
  - "continue BettaMt run"
---

# bettamt-run — BettaMt end-to-end executor (manual-first)

This skill is a thin router plus the explicit execution recipes. It never edits
`params.json` and never invents arguments: everything it runs comes from
`params.json` written by `bettamt-preflight`.

## 0. Resolve skill root, runner and run directory

```bash
: "${BETA_MT_HOME:=$(cd "$SKILL_DIR/../.." && pwd)}"    # skill root = repo root
NF_HOME="$BETA_MT_HOME/runners/nextflow-runner"
test -f "$NF_HOME/main.nf" || { echo "Cannot locate the BettaMt runner tree. Set BETA_MT_HOME."; exit 1; }

# Run dir: wherever outputs should land. Default: $BETA_MT_HOME/runs/<sample>.
: "${RUN_DIR:=$BETA_MT_HOME/runs/default}"
mkdir -p "$RUN_DIR"
```

## 1. Determine the user's stage (filesystem first, then ask)

```bash
# A polished mitogenome exists → run succeeded, time for QC
: "${RUN_DIR:=$PWD}"
if ls "$RUN_DIR"/*_mitogenome.fasta >/dev/null 2>&1; then
    STAGE="qc"
elif [ -f "$RUN_DIR/annotation-report.md" ]; then
    STAGE="done"
elif [ -f "$RUN_DIR"/annotated_genes.fasta ]; then
    STAGE="annotate-qc"
elif [ -f "$RUN_DIR/diagnosis.md" ]; then
    STAGE="review-diagnosis"
elif [ -f "$RUN_DIR/params.json" ]; then
    STAGE="run"
else
    STAGE="preflight"
fi
echo "Stage: $STAGE"
```

If ambiguous, ask one question: *new run, continuing with params.json, failed
run, or success needing QC?*

## 2. Route by stage

| Stage | Action |
|---|---|
| `preflight` | invoke `preflight/bettamt-preflight` (produces `params.json` + `execution_mode`), return here |
| `run` | execute §3 below — recipes (manual, default) or NF runner |
| `debug` | invoke `debug/bettamt-debug` (writes `diagnosis.md`), return |
| `qc` | invoke `qc/bettamt-qc` (writes `report.md`) |
| `annotate-qc` | run `qc/bettamt-qc` §6 annotation checks (writes `annotation-report.md`) |
| `review-diagnosis` | walk the user through Option A in `$RUN_DIR/diagnosis.md` |

**Do not skip preflight.** Even if the user supplies params manually, run
preflight first to validate choices and produce the audit trail.

### a) Execution-mode fork (SP1a — no re-asking)

`params.json:execution_mode` was decided once at preflight. Honour it silently:

| Value | Runner | Source of truth |
|---|---|---|
| `"manual"` (default) | `runners/bash-recipes/*.sh` | the recipe scripts (mirror the NF modules, enforced by `test_smoke.sh`) |
| `"nextflow"` | `nextflow run runners/nextflow-runner/main.nf ...` | the NF module `script:` blocks |

## 3. Execute

### 3a. manual mode — bash recipes (default)

Derive arguments **only** from `params.json` (platform, reads, ref_mito, size,
taxon, rounds, sample_id). Confirm with the user before spending compute:

> "This will run the <ont|illumina> path to a polished mitogenome in `$RUN_DIR`. Proceed? (yes / no / dry-run first)"

```bash
# ONT (reads = single FASTQ)
bash "$BETA_MT_HOME/runners/bash-recipes/run_ont.sh" \
    --reads "<params.reads>" --ref-mito "<params.ref_mito>" \
    --size "<params.size>" --taxon "<params.taxon>" \
    --sample-id "<params.sample_id>" --outdir "$RUN_DIR" --threads "$THREADS"

# Illumina (reads_r2 also required)
bash "$BETA_MT_HOME/runners/bash-recipes/run_illumina.sh" \
    --reads "<params.reads>" --reads-r2 "<params.reads_r2>" \
    --ref-mito "<params.ref_mito>" \
    --rounds "<params.rounds>" --taxon "<params.taxon>" \
    --sample-id "<params.sample_id>" --outdir "$RUN_DIR" --threads "$THREADS"
```

- `--threads` is the only resource knob (preflight picks `min(ncpu, 8)` for a
  workstation unless the machine is small; GetOrganelle/SPAdes and Flye are the
  memory-hungry stages, but a mitogenome fits comfortably in 16 GB).
- The recipes print `[stage k/N] NAME` headers. On failure the last header names
  the failing stage — carry it into `debug/bettamt-debug` as the first evidence.
- Re-run a stage in isolation by copying its command block from the recipe into
  `$RUN_DIR/stages/` — do not ad-hoc-edit the recipe.
- Optional annotation (skipped silently unless `params.ref_gff` or
  `params.ref_gb` is set): after the polished FASTA exists,

```bash
bash "$BETA_MT_HOME/runners/bash-recipes/annotate_genes.sh" \
    --target "$RUN_DIR/<sample>_mitogenome.fasta" \
    --ref-gff "<params.ref_gff>" --ref-fasta "<params.ref_gff_fasta>"   # or --ref-gb "<params.ref_gb>"
```

### 3b. nextflow mode — opt-in runner

Check for the preflight-written user override, then run from the runner tree:

```bash
PROFILE_LOWER=$(echo "${BETA_PROFILE:-local}" | tr '[:upper:]' '[:lower:]')
USER_CONFIG="$RUN_DIR/nextflow.config.${PROFILE_LOWER}.local"
[ -f "$USER_CONFIG" ] && CONFIG_FLAG="-c $USER_CONFIG" || CONFIG_FLAG=""

cd "$BETA_MT_HOME"
nextflow run runners/nextflow-runner/main.nf -profile "$PROFILE_LOWER" $CONFIG_FLAG \
    -params-file "$RUN_DIR/params.json"          # add -dry-run on first use; add -resume to retry
```

- **Never skip the user confirmation** ("OK to proceed? yes/no/dry-run").
- Dry-run first for new users (shows the DAG + resolved per-process resources).
- `-resume` is a separate, distinct command for retries; if the user override
  changed between runs, `-resume` retries the failed task under the new limits.
- The `local` profile inherits `cpus=1 / memory='4 GB'` from `conf/base.config`;
  without an override `GETORGANELLE_ASM` (SPAdes) can OOM — that is why preflight
  writes `nextflow.config.<env>.local`.
- Never run two NF runs in the same `$RUN_DIR` concurrently — `-resume` cannot
  disambiguate them.

## 4. After the run

| Outcome | Next |
|---|---|
| exit 0, `<sample>_mitogenome.fasta` exists | `qc/bettamt-qc` |
| annotation outputs exist (annotated_genes.\*) | `qc/bettamt-qc` §6 — assembly QC first, annotation QC second |
| exit non-zero, no diagnosis yet | `debug/bettamt-debug` |
| `diagnosis.md` exists | walk Option A with the user |

Hand off explicitly — name the next sub-skill and its artifact, don't make the
user discover it.

## 5. Important rules

- **Never** execute without explicit user confirmation.
- **Never** edit `params.json` after preflight — re-run preflight instead.
- **Never** change a recipe or a module one-sidedly: `runners/bash-recipes/test_smoke.sh` enforces the mirror contract and the build gate fails on drift.
- **Always** surface the final artifact path and the next phase at the end.
- **If `$BETA_MT_HOME` doesn't resolve, ask** — a wrong home is wasted compute.

## 6. Output contract — no artifact

This skill produces no file of its own: the runners write `$RUN_DIR`, and all
review artifacts (`report.md`, `annotation-report.md`, `diagnosis.md`) belong to
their phases.

## 7. Common follow-ups

| User says | What to do |
|---|---|
| "Looks good" after QC | summarize: runner mode, polished FASTA path, key QC verdicts |
| "Can I publish this?" | point to `report.md`, `annotation-report.md`, `params.rationale.md` as supplementary material |
| "I want to run more samples" | new `$RUN_DIR` per sample (`$BETA_MT_HOME/runs/<sample>`), same params schema |
| "How do I cite this?" | `cat "$BETA_MT_HOME/README.md"` — citations block at the bottom |
| "It failed at stage 3 (GETORGANELLE_FILTER)" | invoke `debug/bettamt-debug` with that stage name + stderr |