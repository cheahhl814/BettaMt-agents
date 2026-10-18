# BettaMt bash-recipes — the default (manual) runner

Single-machine, single-sample execution path for the BettaMt mitogenome workflow.
These scripts are the **source of truth**: they mirror the Nextflow runner
(`../nextflow-runner/`) command-for-command, and `test_smoke.sh` enforces the
lock-step. The Nextflow runner is the opt-in convenience layer for cluster /
resume-heavy use.

## Scripts

| Script | Mirrors | What it does |
|---|---|---|
| `run_ont.sh` | `subworkflows/ont_mitogenome.nf` | BAIT_MITO → FLYE → GETORGANELLE_FILTER → FIRSTGENE → CIRCLATOR → RACON |
| `run_illumina.sh` | `subworkflows/illumina_mitogenome.nf` | QC_SHORT → GETORGANELLE_ASM → FIRSTGENE → CIRCLATOR → POLCA |
| `annotate_genes.sh` | `modules/local/annotate_genes.nf` | optional post-polish BLAST gene annotation (`bin/blast_extract_cds.py`) |
| `test_smoke.sh` | all | structural lock-step checks (see below) |

## Usage

```bash
# ONT (defaults: --size 16k --taxon vertebrate --threads 4)
runners/bash-recipes/run_ont.sh \
    --reads sample.fastq.gz \
    --ref-mito reference_mitogenome.fasta \
    --sample-id sample --threads 8 --outdir ./runs/sample

# Illumina paired-end (defaults: --rounds 15 --taxon vertebrate)
runners/bash-recipes/run_illumina.sh \
    --reads sample_R1.fastq.gz --reads-r2 sample_R2.fastq.gz \
    --ref-mito reference_mitogenome.fasta \
    --sample-id sample --threads 8

# Optional annotation (GFF3 route, tRNAscan on by default)
runners/bash-recipes/annotate_genes.sh --target runs/sample/sample_mitogenome.fasta \
    --ref-gff ref.gff3 --ref-fasta ref.fasta
```

Both run scripts print per-stage headers (`[stage k/N] NAME`) and finish with a
polished `$outdir/<sample>_mitogenome.fasta`, ready for `qc/bettamt-qc`.

## Environment pinning

All tool invocations go through
`pixi run --manifest-path <skill-root>/runners/nextflow-runner/pixi.toml ...` —
the same pinned env the Nextflow runner uses (one manifest per skill, per the
pixi-skill convention; never per-sub-skill manifests).

Threads default to 4; override with `--threads N`. The ONT/illumina scripts are
single-process-per-stage, so thread count is the only resource knob needed on a
workstation — the Nextflow `conf/*.config` resource blocks apply only to the
NF runner. A mitogenome run fits comfortably in 16 GB RAM.

## Lock-step contract

`test_smoke.sh` greps the module `script:` blocks from `../nextflow-runner/modules/local/*.nf`
and asserts each tool invocation command appears in the mirroring recipe (same
tool name + subcommand, tolerating Groovy interpolation differences). If you
change a recipe or a module, run the smoke test — a one-sided change MUST fail
the check (that is what caught the `.py`-suffix and `_corrected.fasta` drifts on
2026-10-18: see `docs-corpus/README.md` battle notes).

## Known env-name deviations from the NF modules (fixed in both)

- `get_organelle_config.py` / `get_organelle_from_assembly.py` — the pinned env
  ships `.py`-suffixed entry points only; `getorganelle_filter.nf` was patched
  accordingly (fixed together with this recipe in the same commit).
- `pypolca` output `<prefix>_corrected.fasta` (not `<prefix>.fasta`) — patched
  in `polca.nf` + `run_illumina.sh`.
- `fastp` was missing from the runner's pixi.toml deps until this commit
  (QC_SHORT silently resolved a host-global fastp otherwise).