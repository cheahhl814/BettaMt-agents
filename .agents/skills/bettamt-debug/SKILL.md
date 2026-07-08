---
name: bettamt-debug
description: Interpret a failed BettaMt Nextflow run, including OOMs caused by profile-vs-machine mismatch. Use when nextflow log shows ERR, exit status 1, exit 137 (OOM-kill), or a process was OOM-killed. Triggers: "bettamt failed", "BettaMt error", "why did BettaMt fail", "nextflow failure", "bettamt debug", "GetOrganelle no valid assembly graph", "Flye failed", "Racon OOM", "out of memory", "exit 137", "spades OOM".
---

# bettamt-debug — BettaMt failure interpreter

This skill is **reactive** — invoked after a `nextflow run` has produced an error. It does not run the pipeline.

## 0. Locate the pipeline and the run

```bash
: "${BETA_MT_HOME:=$(cd "$SKILL_DIR/../../../BettaMt" 2>/dev/null && pwd)}"
test -f "$BETA_MT_HOME/main.nf" || { echo "Cannot locate BettaMt. Set BETA_MT_HOME."; exit 1; }

# Run dir is wherever the user launched nextflow (usually $BETA_MT_HOME or ./results)
RUN_DIR="${RUN_DIR:-$BETA_MT_HOME}"
test -f "$RUN_DIR/.nextflow.log" || test -f "$RUN_DIR/nextflow.log" || { echo "No nextflow.log found in $RUN_DIR"; exit 1; }
```

Identify the run id from the latest trace:

```bash
RUN_ID=$(nextflow log -q "$RUN_DIR" | tail -1)
```

## 1. Identify the failing process(es)

```bash
nextflow log -f "process,name,status,exit,hash,duration" "$RUN_DIR" \
  | awk -F'\t' '$3 != "COMPLETED" { print }'
```

You'll get a list of failing tasks with their work-dir hashes. For each failure:

```bash
# Get the work directory
nextflow log -f "process,name,exit,hash,workdir" "$RUN_DIR" | grep -F "$HASH"
```

Then read:

- `$WORKDIR/.command.err` — the actual stderr from the failing command
- `$WORKDIR/.command.log` — the full stdout
- `$WORKDIR/.command.sh` — the exact command that was run (useful for reproducing)
- `$WORKDIR/.exitcode` — the exit code

## 1.5 OOM diagnosis: profile vs reality

Before falling into the per-process signature library, check whether the OOM is a **config-vs-machine** problem. This is the single most common cause of OOM in BettaMt runs, and the fix is mechanical (not diagnostic):

```bash
# 1. What did the profile / override ask for?
HASH=<failing hash>
nextflow log -f "process,name,cpus,memory,realtime,hash" "$RUN_DIR" | grep -F "$HASH"

# 2. How much memory does the host have?
free -g | head -3

# 3. Was a user-side override applied?
USER_CONFIG="$RUN_DIR/nextflow.config.${BETA_PROFILE:-local}.local"
[ -f "$USER_CONFIG" ] && echo "User override in use: $USER_CONFIG" \
                      || echo "No user override — using $BETA_MT_HOME/conf/{base,local,slurm}.config defaults"
echo ""
echo "BettaMt base.config defaults: cpus=1, memory='4 GB' for ALL processes"
echo "The local profile inherits these. Without a user override, even on a 64 GB"
echo "laptop, every BettaMt process gets 1 cpu / 4 GB — which OOMs on GetOrganelle_ASM."
```

**If the process requested more memory than the host has, the answer is the override file, not a tool fix.** Generate or update `nextflow.config.<env>.local` (see `bettamt-preflight` §2.5) and `-resume`.

**Pattern**: `oom-kill` exit code is `137` on Linux; on SLURM look for `slurmstepd: ... OOM` in `slurm-<jobid>.out`. The Nextflow task trace file `<WORKDIR>/.command.trace` shows the requested vs used memory:

```bash
cat <WORKDIR>/.command.trace
# Look at:  memory, peak_memory, vmem, rss
```

If `vmem` > host RAM, it's a config problem. If `vmem` < host RAM but the process still died with OOM, it's a real tool bug or a leak — fall into §2's signature library.

**Special-case for BettaMt's `local` profile**: because `conf/local.config` is just `executor = 'local'`, every process inherits the `base.config` defaults of 1 cpu / 4 GB. So **a `local` user without an override will OOM on `GETORGANELLE_ASM` (SPAdes) and `FLYE` almost every time**, regardless of host RAM. The override isn't optional for local runs on heavy processes.

## 2. Signature library (extend as you encounter new failures)

Match the failing command's stderr / exit code against these known patterns. **Always** read the actual error before concluding — never pattern-match blindly.

### GetOrganelle (short-read)

| Signature in `.command.err` | Likely cause | Suggested fix |
|---|---|---|
| `No valid assembly graph found!` | SPAdes graph failed disentanglement; usually low coverage or very high coverage | If coverage was low: more sequencing. If very high: re-run with `--reduce-reads-for-coverage 2000` (override in `modules/local/getorganelle_asm.nf`) |
| `Assembling exited halfway` | SPAdes OOM or crashed | Bump `withName:'GETORGANELLE_ASM' { memory = '128 GB' }` in your user-side override (`nextflow.config.<env>.local`, see `bettamt-preflight` §2.5). Do **not** edit `$BETA_MT_HOME/conf/slurm.config`. |
| `No paired reads found?!` | R1/R2 are swapped, or files are corrupted | Swap R1↔R2; verify with `seqkit head -n 4 R1.fq.gz R2.fq.gz` and confirm read names match |
| `mean error rate = 0.4+` (very early in log) | Wrong encoding (color-space / old Phred+64) | Re-export reads as Illumina 1.8+; check `head -c 100 R1.fq.gz` |
| `seed bowtie2 index existed! ... 0 seed reads` | `--ref_mito` not a mitogenome or wrong taxon | Verify with `tRNAscan-SE` and `BLASTn` against NCBI organelle DB |
| `Hit the round N and stop` followed by `No valid` | `--rounds` too low for coverage / genome complexity | Increase `--rounds` (preflight default 15, but try 30) |

### Flye (long-read)

| Signature | Likely cause | Suggested fix |
|---|---|---|
| `ERROR: Inconsistent edge length` | Mix of read chemistries | Use only one chemistry; re-basecall |
| `killed` (after long runtime) | OOM | Bump `withName:'FLYE' { memory = '128 GB' }` in your user-side override; reduce `--genome-size` to actual mitogenome size |
| `assembled contigs length: 0` | No reads passed baiting | Verify `bait_mito` output; lower `bbduk hdist` to 2; relax seqkit length filter |
| `Oops, something went wrong` (no traceback) | Corrupted FASTQ | Re-download / re-basecall; check FASTQ with `seqkit stats` |

### Racon (ONT polishing)

| Signature | Likely cause | Suggested fix |
|---|---|---|
| `Segmentation fault` | Insufficient memory or huge alignments | Bump `withName:'RACON' { memory = '...' }` in your user-side override; subset reads with `seqkit sample -n 50000` first |
| `Cannot find minimap2 in $PATH` | Env not activated | Verify `pixi run --manifest-path $BETA_MT_HOME/pixi.toml which minimap2` |

### POLCA / pypolca (Illumina polishing)

| Signature | Likely cause | Suggested fix |
|---|---|---|
| `freebayes: command not found` | pixi env didn't install `freebayes` (transitive dep of `pypolca`) | `pixi add freebayes` in `$BETA_MT_HOME/pixi.toml` |
| `BamFile ... not indexed` | `samtools` < 1.18 conflict | Confirm `pixi.toml` pins `samtools >= 1.18` (we already do) |
| `Assembly has 0 contigs` | Upstream CIRCLATOR failed | Check `results/circlator/` — usually means tRNA-Phe wasn't found |

### Shared / annotation (FIRSTGENE, CIRCLATOR)

| Signature | Likely cause | Suggested fix |
|---|---|---|
| `No tRNA-Phe found` | `--taxon` model wrong for this lineage (e.g. `mammal` on a fish mitogenome) | Switch `--taxon vertebrate` (or `anonym` as a fallback) |
| `trnF.bed: empty` | tRNAscan found no Phe, or contig is not a mitogenome | Re-run `tRNAscan-SE` on the raw contig; verify contig length is in 14–20 kb window |
| `circlator fixstart: cannot find genes_fa` | Upstream FIRSTGENE failed | Re-run with `-resume` after fixing FIRSTGENE |
| `circlator clean: contigs shorter than --min_contig_length` | Mitogenome is fragmented | Lower `--min_contig_length 12000` to `8000` in `modules/local/getorganelle_filter.nf` (only if you trust the input) |

### Gene annotation (ANNOTATE_GENES / blast_extract_cds.py)

Only runs if `--ref_gff` or `--ref_gb` was supplied.

| Signature in `.command.err` | Likely cause | Suggested fix |
|---|---|---|
| `blastn: command not found` | pixi env missing `blast` | `pixi add blast` in `$BETA_MT_HOME/pixi.toml` (currently `>=2.14`) |
| `tRNAscan-SE: command not found` | tRNAscan-SE not on PATH; the script's auto-detect in `TRNASCAN_SEARCH_PATHS` missed the pixi env | Pass `--trnascan` explicitly to the run, or set the pixi env as the active shell |
| `WARNING: reference missing: ...` | Reference GFF/GB doesn't have all 13 PCGs (common for incomplete annotations) | Use a more complete reference; warnings are non-fatal and the present genes still annotate correctly |
| `No D-loop detected (no gap >= 100 bp)` | Mitogenome is very compact; the D-loop might be < 100 bp or the assembly didn't include it | Verify the polished contig actually contains the D-loop region (length should be 14–20 kb); if shorter, the assembly is truncated |
| `FileNotFoundError: --ref-fasta` | Used `--ref-gff` without `--ref-fasta` and the GFF3 has no `##FASTA` section | Add `--ref_gff_fasta` to params.json or switch to `--ref_gb` (self-contained) |
| `No BLAST hit` for most genes | Reference is too divergent (different order/family) | Use a closer reference; consider dropping annotation entirely if no suitable reference exists |
| `BLAST hit coordinates overlap / duplicated gene` | The polished mitogenome has genuine gene duplications (NUMT or recent duplication) | Inspect with IGV; this is biological, not a pipeline bug |
| `KeyError: 'D-loop'` or similar from blast_extract_cds.py | Bug in gene name normalization; the `GENE_ALIASES` dict should be extended | Open a bug report with the offending reference file |

### Resource / scheduler

| Signature in `nextflow.log` | Likely cause | Suggested fix |
|---|---|---|
| `Process requirement exceeds available CPUs` | SLURM `--cpus` bigger than partition limit | Reduce `cpus` for that process in your user-side override (`nextflow.config.<env>.local`). Do **not** edit `$BETA_MT_HOME/conf/slurm.config`. |
| `slurmstepd: Exceeded job memory limit` | OOM-killed by SLURM | Bump `memory` directive for that process in your user-side override |
| `Exit status 137` | SIGKILL = OOM | **First** check §1.5 — this is almost always a config-vs-machine mismatch in BettaMt. Bump memory in the override (if headroom) or push the run to a larger node |
| `Exit status 143` | SIGTERM = timeout or manual cancel | Re-run with `-resume` |
| Hardcoded `--partition=cpu1.medm` rejected by SLURM | Partition name doesn't exist on this cluster | Override in `nextflow.config.slurm.local`; the user's cluster has different partition names. See `bettamt-preflight` §2.5d. |

## 3. Output contract — write `diagnosis.md`

Always write to `$RUN_DIR/diagnosis.md`. Template:

```markdown
# BettaMt failure diagnosis

Run:     $RUN_ID
When:    <ISO8601>
Profile: slurm|local
Pipeline: $BETA_MT_HOME

## Failed processes (in order of execution)

### 1. PROCESS_NAME (exit N)
- **Work dir**: $WORKDIR
- **Error excerpt**:
  \`\`\`
  <3-10 lines from .command.err>
  \`\`\`
- **Likely cause** (ranked, with confidence):
  1. **<hypothesis>** — <one-line explanation>. *(high confidence — matches signature '...' in skill library)*
  2. **<hypothesis>** — <one-line explanation>. *(medium)*
  3. **<hypothesis>** — <one-line explanation>. *(low)*

## Recommended next action (pick one)

### Option A — cheapest fix, try first
\`\`\`bash
<exact command the user can paste>
\`\`\`

### Option B — if A didn't help
\`\`\`bash
<exact command>
\`\`\`

## Diagnostic artifacts to inspect manually
- $WORKDIR/.command.err
- $WORKDIR/.command.log
- $WORKDIR/.command.sh
- $RUN_DIR/.nextflow.log

## Skill library extension
- If you saw a *new* error pattern, add it to the table in
  `bettamt-debug/SKILL.md` so future runs are diagnosed faster.
```

Always produce **multiple ranked hypotheses** (top 3) — pipelines fail in many ways, and picking the first match is how you get stuck in local minima. Use the confidence labels explicitly.

## 4. What NOT to do

- Do **not** modify `main.nf`, modules, or config without the user's explicit go-ahead
- Do **not** blanket-edit `$BETA_MT_HOME/conf/slurm.config` to "try bigger numbers" — diagnose the actual OOM first (§1.5), then update the user-side `nextflow.config.<env>.local` override
- Do **not** use `withLabel:` in a BettaMt override — BettaMt targets processes by `withName:` (e.g. `'BAIT_MITO'`, `'GETORGANELLE_ASM'`); a `withLabel:` block will silently do nothing
- Do **not** clear the `work/` directory without checking if `-resume` would have been enough
- Do **not** pattern-match a hypothesis without quoting the actual error line
- Do **not** claim a single hypothesis is the only one — always give the top 3

## 5. Resume vs re-run

If the failure is **deterministic** (config bug, wrong params, missing file) → fix the cause, then `nextflow run -resume` (picks up cached steps).

If the failure is **stochastic** (transient OOM, network blip, scheduler hiccup) → `nextflow run -resume` is often enough.

If the failure is **upstream** (e.g. `bait_mito` produced nothing) → `-resume` won't help; fix upstream and re-run from scratch.

State which case applies in `diagnosis.md`.
