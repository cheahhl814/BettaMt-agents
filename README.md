# BettaMt-agents — agentic mitogenome assembly, annotation & QC

> **Skill bundle v1.0.0** — rebuilt to the skill-forge convention (root router +
> phase sub-skills + dual runners) with **manual bash recipes as the default
> execution path** and the Nextflow pipeline demoted to an opt-in cluster runner.
> Predecessor layout: `.agents/skills/*` (legacy BettaMt-pattern skills, now
> folded into `preflight/ run/ qc/ debug/` under the same commit history).

![agents](https://img.shields.io/badge/agent-skill--forge-blue) ![version](https://img.shields.io/badge/version-1.0.0-blue) ![platform](https://img.shields.io/badge/platform-linux--64-green)

## Contents

- [Pipeline overview](#-pipeline-overview)
- [Installation](#-installation)
- [Usage](#-usage)
- [Tools](#-tools)
- [Update check](#-update-check)
- [Repository layout](#-repository-layout)
- [Hard guarantees](#-hard-guarantees)
- [Provenance](#provenance)

## 🔬 Pipeline overview

Assembles a vertebrate **mitochondrial genome** from raw reads plus a reference
mitogenome, polishes it, optionally transfers gene annotations, and QC-reports
the result — with the agent guided phase-by-phase
(preflight → run → qc → debug) and every non-obvious decision gated behind an
Evidence → Recommend → Options ask-user stop point.

| Path | Reads | Stages |
|---|---|---|
| **ONT** (long reads) | single FASTQ | bait & remap (seqkit → bbduk → minimap2 → samtools) → Flye `--meta` assembly → GetOrganelle graph extraction + circlator clean → trnF extraction (tRNAscan-SE) → circlator fixstart → racon polish |
| **Illumina** (paired reads) | R1 + R2 | sequali + fastp trim → GetOrganelle seed-guided assembly → trnF/circlator fixstart → pypolca polish |
| optional | — | BLAST-based gene annotation (`bin/blast_extract_cds.py`, GFF3 or GenBank reference) on the polished mitogenome |

## 🚀 Installation

```bash
git clone https://github.com/cheahhl814/BettaMt-agents
cd BettaMt-agents
pixi install --manifest-path runners/nextflow-runner/pixi.toml   # pre-pull the pinned tool env
pixi install --manifest-path runners/nextflow-runner/envs/nextflow/pixi.toml   # only for the NF-mode smoke
# agents: symlink-deploy the SKILL ROOT so the deployed copy is the live tree:
ln -sfn "$PWD" ~/.agents/skills/bettamt
```

## 💡 Usage

```bash
# 1. Preflight: validate inputs, pick platform/taxon/size/rounds, write params.json
#    (agent: invoke preflight/bettamt-preflight/)
# 2. Run — manual mode is the default (single machine, per-stage visibility):
bash runners/bash-recipes/run_ont.sh \
    --reads my_sample.fastq.gz --ref-mito reference_mito.fasta \
    --sample-id my_sample --threads 8 --outdir ./runs/my_sample
# or Illumina:
bash runners/bash-recipes/run_illumina.sh \
    --reads my_R1.fastq.gz --reads-r2 my_R2.fastq.gz \
    --ref-mito reference_mito.fasta --rounds 15 --sample-id my_sample
# 3. QC the polished mitogenome (agent: invoke qc/bettamt-qc/) -> report.md
# 4. Optional annotation after a PASS:
bash runners/bash-recipes/annotate_genes.sh --target runs/my_sample/my_sample_mitogenome.fasta \
    --ref-gff ref.gff3 --ref-fasta ref.fasta
```

Nextflow mode (opt-in — cluster / resume / trace reports):

```bash
nextflow run runners/nextflow-runner/main.nf -profile slurm \
    -c runs/<sample>/nextflow.config.slurm.local \
    -params-file runs/<sample>/params.json
```

`bettamt-preflight` writes the user-side `nextflow.config.<env>.local` override
(never edit `conf/*.config` in place); `run/bettamt-run` documents the
`-resume` and `-dry-run` conventions.

## 🧰 Tools

Every tool invocation goes through
`pixi run --manifest-path runners/nextflow-runner/pixi.toml` — **one pinned
manifest shared by both runners** (per the pixi-skill convention: never fork
the env per sub-skill). Pinned: seqkit, bbmap (bbduk.sh), minimap2, samtools,
flye, getorganelle (`.py` entry points), circlator, racon, bedtools,
trnascan-se, fastp, sequali, pypolca, blast. Verbatim CLI references for every
pinned version live in `docs-corpus/<tool>/cli-reference.md`.

An **external agent-side dependency** (optional, degrades gracefully):
[read-qc-trimming](https://github.com/cheahhl814/read-qc-trimming) — preflight
can route reads through its QC→trim→QC sandwich before BettaMt consumes them;
otherwise the pinned in-run `QC_SHORT` stage (sequali + fastp) applies.

## 🔄 Update check

```bash
pixi run update-check
```

`bin/skill-update-check.py` diffs the deployed copy's git SHA + SKILL.md version
against `origin` (`git fetch`, no API calls — works air-gapped). Verdicts:
`UP-TO-DATE` / `LOCAL-AHEAD` / `BEHIND-BY-N` (exit 1 → re-deploy) / `OFFLINE` /
`NO-ORIGIN` (exit 2, informational). Re-deploy convention: symlink the repo
root (see Installation) and verify with `diff -rq`.

## 📁 Repository layout

```
BettaMt-agents/                 # $BETTAMT_HOME = skill root = repo root
├── SKILL.md                    # agent entry point (router)
├── preflight/bettamt-preflight/  # params.json + params.rationale.md + execution_mode
├── run/bettamt-run/              # execution-mode fork (manual default | nextflow opt-in)
├── qc/bettamt-qc/                # assembly QC + annotation QC (§6)
├── debug/bettamt-debug/          # signature-library failure interpreter
├── runners/
│   ├── bash-recipes/             # DEFAULT runner: run_ont.sh, run_illumina.sh, annotate_genes.sh
│   │                             #   + test_smoke.sh (recipe↔module lock-step contract)
│   └── nextflow-runner/          # opt-in runner: DSL2 pipeline + conf/{local,slurm}
│       │                         #   + pixi.toml/pixi.lock (the pinned env) + bin/blast_extract_cds.py
│       └── envs/nextflow/        # nextflow runtime for the smoke harness only
├── docs-corpus/                  # pinned --help snapshots for every tool + battle notes
├── pixi.toml                     # skill-level tasks (update-check, smoke)
└── bin/skill-update-check.py     # deployed-copy vs upstream SHA diff
```

## 🔒 Hard guarantees

- **Manual recipes are the source of truth.** They mirror the Nextflow module
  scripts command-for-command; `runners/bash-recipes/test_smoke.sh` enforces the
  lock-step, so a one-sided change fails the build.
- **One pinned env for both runners.** Tool versions are reproducible via
  `pixi.lock`; no host-global fallback (this exact class of bug — host-global
  `fastp` silently reading 0 reads — was found and fixed at rebuild time).
- **No invented flags or model names.** Agents read `docs-corpus/` first.
- **No silent destructive actions.** Execution always confirms with the user;
  ask-user stop points are encoded per phase (SP tables).
- **No GPU anywhere.** The whole workflow is CPU tools; nothing recommends GPU.
- **Run artifacts never ship.** `runs/`, `params.rationale.md`,
  `nextflow.config.*.local`, `report.md`, `diagnosis.md` etc. are gitignored.

## Provenance

- **BettaMt** Nextflow pipeline + legacy agent skills — this repo's own history
  (commits `3eafa7f`…`f8fd3e4`), rebuilt under the same git history.
- **skill-forge** convention and `bin/skill-lint.py` —
  github.com/cheahhl814/skill-forge; layout mirrors
  `bacterial-genome-analysis` v5.x and the nanobc v2.x dual-runner pattern.
- The mitogenome toolchain owes everything to Flye, GetOrganelle, Circlator,
  racon, pypolca, tRNAscan-SE, bbmap, minimap2, samtools, seqkit, bedtools,
  fastp, sequali and BLAST+ — see `docs-corpus/<tool>/cli-reference.md` for the
  pinned versions this bundle uses. Citation: `CITATION.cff`.