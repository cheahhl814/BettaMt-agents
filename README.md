# BettaMt

[![Version](https://img.shields.io/badge/version-0.2.0-blue)](#-installation)
[![Type](https://img.shields.io/badge/type-agent%20skill-blueviolet)](#-installation)
[![Built with](https://img.shields.io/badge/built%20with-bioinfo--skill--creator-orange)](https://github.com/cheahhl814/bioinfo-skill-creator)
[![Shape](https://img.shields.io/badge/shape-multi--tool%20pipeline-8A2BE2)](#-pipeline-overview)

A [Nextflow](https://www.nextflow.io/) DSL2 pipeline that assembles and polishes mitochondrial genomes from either Oxford Nanopore long reads or Illumina paired-end short reads, with every tool resolved through a reproducible [pixi](https://pixi.sh) environment. It ships with a companion layer of [Agent Skills](https://agentskills.io/specification) that drive the pipeline end-to-end — parameter preflight, run orchestration, failure interpretation, and post-assembly / post-annotation quality control — following the AiX-BIO skill convention.

**Repository**: https://github.com/cheahhl814/BettaMt-agents

> [!NOTE]
> Current version: **v0.2.0** (updated 2026-09-17).

## Contents

- [Installation](#-installation)
- [Usage](#-usage)
- [Pipeline overview](#-pipeline-overview)
- [Tools](#-tools)
- [Update check](#-update-check)
- [Repository layout](#-repository-layout)
- [Hard guarantees](#-hard-guarantees)
- [Provenance](#-provenance)

## 🚀 Installation

This is an **agent-driven pipeline**, not a user-facing library. The recommended install path is to let your AI agent import it.

**Option A — give your agent this prompt (recommended):**

```text
Install the BettaMt mitogenome pipeline from
https://github.com/cheahhl814/BettaMt-agents —
clone it, run `pixi install` inside the BettaMt/ directory, make sure
`nextflow` is on your PATH, and read the SKILL.md files under
.agents/skills/ to understand the phase chain and each skill's
Input/Output contract. Confirm when the environment is ready.
```

**Option B — manual install:**

```bash
git clone https://github.com/cheahhl814/BettaMt-agents.git
cd BettaMt-agents/BettaMt
pixi install               # resolves every pinned tool from pixi.toml
pixi global add nextflow   # or use any system Nextflow ≥ 23.10
```

> [!TIP]
> `pixi install` resolves every pinned tool from `pixi.toml` (channels: conda-forge, bioconda) into an isolated `.pixi/` environment — no system-wide installs, no version conflicts with other skills.

The skills under `.agents/skills/` follow the [Agent Skills standard](https://agentskills.io/specification) and are auto-discovered project-locally by agents that scan `.agents/skills/` (pi, OpenCode, Codex). For agents that don't scan that path (e.g. Claude Code), add a symlink in the project root: `ln -s .agents/skills .claude/skills`, or copy the folder into your agent's skills directory.

## 💡 Usage

The workflow is designed to be driven by an AI agent: the agent reads the orchestrator skill (`bettamt-run`), detects the current phase from the filesystem state of your run directory, and routes to the right sub-skill.

### Natural-language prompts that trigger the skill

```text
run BettaMt / assemble my mitogenome
set up BettaMt / what params for BettaMt
BettaMt failed / why did BettaMt fail / exit 137
check my mitogenome / did the assembly work
check my gene annotation / did the BLAST annotation work
```

### Manual phase execution

If you prefer to drive the phases yourself (or want to re-run a single phase):

```bash
# Set your run directory (all handoff artifacts land here)
RUN_DIR=/path/to/run-dir

# Phase 1 — preflight (skill: bettamt-preflight): produce params.json
#   + params.rationale.md from your FASTQ + reference mitogenome.

# Phase 2 — run the pipeline (skill: bettamt-run), after preflight ≥ GO:
cd BettaMt

# ONT (long reads)
nextflow run main.nf -profile local \
  --platform ont \
  --reads    data/reads.fastq \
  --ref_mito data/reference.fasta \
  --size     16k \
  --taxon    vertebrate

# Illumina (paired-end short reads)
nextflow run main.nf -profile local \
  --platform   illumina \
  --reads      data/sample_R1.fastq.gz \
  --reads_r2   data/sample_R2.fastq.gz \
  --ref_mito   data/reference.fasta \
  --rounds     15 \
  --sample_id  betta_splendens \
  --taxon      vertebrate

# Phase 3 — on failure, interpret stderr via the signature library
#   (skill: bettamt-debug) → diagnosis.md

# Phase 4 — QC the outputs (skills: bettamt-qc, bettamt-annotate-qc)
#   → report.md, annotation-report.md
```

> [!IMPORTANT]
> Each phase has an explicit **Input/Output contract** at the top of its `SKILL.md`. If the upstream artifact is missing (e.g. you run before preflight passed), the skill refuses to proceed and tells you which phase to run first. Phases are gated on purpose — don't skip them.

> [!NOTE]
> The skills resolve the pipeline directory via the `BETA_MT_HOME` environment variable (default: the `BettaMt/` sibling of the skills folder) and write handoff artifacts to `RUN_DIR` (default: next to the pipeline) — the pipeline's own working directory stays clean.

## 🔬 Pipeline overview

The analysis follows a phased evidence chain. Each phase consumes the artifacts of the previous phase and gates progression with a Go/No-Go check.

| # | Phase | Goal | Skill | Artifact produced |
|:--|:------|:-----|:------|:------------------|
| **1** | **Preflight** | Validate FASTQ + reference, choose platform/rounds/taxon, produce Nextflow params + rationale | `bettamt-preflight/` | `params.json`, `params.rationale.md` |
| **2** | **Run** | End-to-end orchestration: detect the user's stage, route, execute `nextflow run` | `bettamt-run/` | Polished mitogenome under `results/` |
| **3** | **Debug** | Interpret a failed run from `.nextflow.log` / `work/` via the failure-signature library | `bettamt-debug/` | `diagnosis.md` |
| **4** | **Assembly QC** | Verdict on completeness, circularity, NUMTs, and publishability | `bettamt-qc/` | `report.md` |
| **5** | **Annotation QC** | Verify the `blast_extract_cds.py` output (`.fasta` / `.gff3` / `.bed`) | `bettamt-annotate-qc/` | `annotation-report.md` |

**Long-read path (`--platform ont`):** `BAIT_MITO` (length filter + map to reference) → `FLYE` assembly → `GETORGANELLE_FILTER` (mito contig extraction) → `FIRSTGENE` (trnF boundary) → `CIRCLATOR fixstart` (circularization + rotation) → `RACON` polish.

**Short-read path (`--platform illumina`):** `QC_SHORT` (QC + adapter/quality trim) → `GETORGANELLE_ASM` (de Bruijn assembly, reference-seeded) → `FIRSTGENE` → `CIRCLATOR fixstart` → `POLCA` polish.

Both paths optionally finish with `ANNOTATE_GENES` (`blast_extract_cds.py`) when `--ref_gff` / `--ref_gb` is supplied, producing `.fasta` / `.gff3` / `.bed` for 13 PCGs, 2 rRNAs, 22 tRNAs, and the D-loop.

> [!TIP]
> Read each skill's `SKILL.md` for the full procedure, its ask-user stop points, and its signature library (stderr pattern → cause → fix).

## 🧰 Tools

All tools are resolved from conda-forge/bioconda via the pinned `BettaMt/pixi.toml`; Nextflow itself is the only separate install.

| Tool | Version | Role |
|:-----|:--------|:-----|
| `nextflow` | ≥ 23.10 | Workflow executor (install separately) |
| `flye` | * | ONT long-read assembly |
| `getorganelle` | ≥ 1.7.7.1, < 2 | Mitochondrial contig extraction (ONT) / de novo assembly (Illumina) |
| `circlator` | * | Circularization, `trnF` fixstart rotation, assembly-graph cleaning |
| `racon` | * | ONT raw-read polishing |
| `pypolca` | ≥ 0.4.0 | Illumina polishing (POLCA) |
| `bwa` | ≥ 0.7.17 | Alignment for polishing |
| `minimap2` | * | Read alignment |
| `bbmap` | * | Read filtering and baiting |
| `samtools` | ≥ 1.18 | Alignment utilities |
| `seqkit` | * | Sequence utilities |
| `blast` | ≥ 2.14 | Gene-annotation transfer (`blast_extract_cds.py`) |
| `trnascan-se` | * | tRNA verification |
| `bedtools` | * | Genomic interval operations |
| `python` | 3.11.6 | Runtime for helper scripts |

Tool invocations are wired into every Nextflow process through pixi (`pixi run --manifest-path ${baseDir}/pixi.toml …`), so flags always match the pinned tool versions.

## 🔄 Update check

This repository does not ship the per-skill self-update check (`bin/skill-update-check.py`) used by skills generated with [bioinfo-skill-creator](https://github.com/cheahhl814/bioinfo-skill-creator): the pipeline and its skills are co-versioned in one repo, so upstream drift is detected directly with git.

```bash
git fetch origin
git rev-list --count HEAD..origin/main   # > 0 → upstream is ahead; git pull to update
```

| State | Meaning |
|:------|:--------|
| `0` | Local checkout matches `origin/main` |
| `N` (> 0) | Upstream is N commits ahead — `git pull` to update |

## 📁 Repository layout

```text
BettaMt-agents/
├── BettaMt/                        # Nextflow DSL2 pipeline (pixi-managed)
│   ├── main.nf
│   ├── bin/
│   │   └── blast_extract_cds.py    # gene-annotation helper
│   ├── conf/ · modules/ · subworkflows/
│   ├── nextflow.config
│   └── pixi.toml                   # pinned tool environment (pixi install)
├── .agents/skills/                 # the agentic layer (one SKILL.md per skill)
│   ├── bettamt-run/                # orchestrator / router (start here)
│   ├── bettamt-preflight/          # params.json + params.rationale.md
│   ├── bettamt-debug/              # failure interpreter
│   ├── bettamt-qc/                 # post-assembly QC
│   └── bettamt-annotate-qc/        # post-annotation QC
├── CITATION.cff                    # GitHub "Cite this repository" metadata
└── .gitignore                      # excludes Nextflow work/, .nextflow/, test data
```

## 🔒 Hard guarantees

- **Filesystem evidence chain** — agent → pipeline communication is `params.json` (with `params.rationale.md` as its rationale); pipeline → agent communication is `.nextflow.log`, `work/`, and `results/`. The pipeline never reads the skill Markdown; the skills never parse `.nf` internals.
- **The agent reasons; the pipeline executes** — new failure modes become new rows in the signature library, not new scripts.
- **Explicit stop points** — ambiguous decisions are surfaced as *Evidence + Recommend + Options*, not auto-picked.
- **Reproducible environments** — every tool is pinned in `pixi.toml` and resolved via pixi.
- **Clean working directory** — handoff artifacts (`params.json`, `diagnosis.md`, `report.md`, `annotation-report.md`) are written next to the pipeline, never inside it; `RUN_DIR` relocates them.

## Provenance

The agentic layer follows the [bioinfo-skill-creator](https://github.com/cheahhl814/bioinfo-skill-creator) skill convention — a phased evidence chain (preflight → run → QC → debug), signature-library debugging, and explicit ask-user stop points. Pattern adopted from:

- **BettaMt-agents** — https://github.com/cheahhl814/BettaMt-agents (origin of the `bettamt-*` pattern)
- **bacterial-genome-analysis** — https://github.com/cheahhl814/bacterial-genome-analysis
- **amr-gene-screening** — https://github.com/cheahhl814/amr-gene-screening

If you use BettaMt in published research, see [`CITATION.cff`](CITATION.cff).

## License

No license file is included yet (set to "LicenseNotSet" in `CITATION.cff`) — add a `LICENSE` (MIT / Apache-2.0 / GPL-3.0) before distributing.