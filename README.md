# BettaMt — Mitogenome assembly pipeline and agentic skills

[![Version](https://img.shields.io/badge/version-0.2.0-blue)](#-installation)
[![Type](https://img.shields.io/badge/type-nextflow%20pipeline%20%2B%20agent%20skills-blueviolet)](#-pipeline-overview)
[![Agent Skills](https://img.shields.io/badge/skills-agent%20skills%20standard-teal)](https://agentskills.io/specification)

A Nextflow DSL2 pipeline for assembling and polishing mitogenomes from long-read (ONT) and short-read (Illumina) sequencing data, paired with five agentic skills that prepare parameters, run the pipeline, QC the assembly, interpret failures, and check annotation output.

**Repository**: https://github.com/cheahhl814/BettaMt-agents

> [!NOTE]
> Current version: **v0.2.0** (updated 2026-09-17).

## Contents

- [Installation](#-installation)
- [Usage](#-usage)
- [Pipeline overview](#-pipeline-overview)
- [Tools](#-tools)
- [Repository layout](#-repository-layout)
- [Hard guarantees](#-hard-guarantees)
- [Provenance](#provenance)

## 🚀 Installation

This repo pairs a pipeline with **agent skills** — Markdown capability packages that an AI agent loads and follows. The recommended install path is to let your AI agent set it up.

**Option A — give your agent this prompt (recommended):**

```text
Install the BettaMt agentic skills from
https://github.com/cheahhl814/BettaMt-agents —
clone the repository into a project directory your agent scans for
`.agents/skills/`, then run `pixi install` inside the `BettaMt/` directory
to resolve the pinned pipeline tools. Read each skill's SKILL.md to
understand its role and its Input/Output contract, and confirm when the
environment is ready.
```

**Option B — manual install:**

```bash
git clone https://github.com/cheahhl814/BettaMt-agents.git
cd BettaMt-agents/BettaMt
pixi install
```

> [!TIP]
> `pixi install` resolves every pinned tool from `pixi.toml` (channels: conda-forge, bioconda) into an isolated `.pixi/` environment — no system-wide installs, no version conflicts with other pipelines.

## 💡 Usage

The five skills are **model-invoked**: the agent reads each skill's `description` frontmatter at session start, and a user prompt matching a trigger phrase causes the agent to load the skill and follow it. Project-local `.agents/skills/` are auto-discovered from your working directory up to the repo root.

### Natural-language prompts that trigger the skills

| Skill | Example prompts | What it produces |
|:------|:----------------|:-----------------|
| `bettamt-preflight` | *"set up BettaMt"*, *"what params for BettaMt"*, *"scale down BettaMt for laptop"* | `params.json` + `params.rationale.md` + machine-specific `nextflow.config` override |
| `bettamt-run` | *"run BettaMt"*, *"assemble my mitogenome"*, *"I have reads and a reference"* | A completed Nextflow run (walks preflight → run → QC/debug) |
| `bettamt-debug` | *"BettaMt failed"*, *"why did BettaMt fail"*, *"exit 137"*, *"Flye failed"* | `diagnosis.md` with top-3 hypotheses |
| `bettamt-qc` | *"check my mitogenome"*, *"is the mitogenome complete"*, *"post-assembly QC"* | `report.md` with pass/warn/fail verdicts |
| `bettamt-annotate-qc` | *"check annotation"*, *"did the BLAST annotation work"*, *"annotate QC"* | `annotation-report.md` |

### Running the pipeline manually

If you prefer to drive the pipeline yourself:

```bash
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
```

`pixi install` is required once on the host so each Nextflow process can resolve its tools via `pixi run --manifest-path ${baseDir}/pixi.toml` from inside its `script:` block — you do **not** wrap the `nextflow` invocation itself in `pixi run`. Put `nextflow` on your `PATH` (e.g. `pixi global add nextflow`).

> [!TIP]
> `--size` is the estimated mitogenome size (`16k` typical for *Betta*); `--taxon` selects the tRNAscan-SE model (`mammal` / `vertebrate`); `--rounds` are GetOrganelle extension rounds (10–15 typical, Illumina only). Full parameter tables: [`BettaMt/README.md`](BettaMt/README.md). Profiles: `local` (default) and `slurm` (explicit per-process resources).

## 🔬 Pipeline overview

The agentic layer follows a phased evidence chain. Each phase consumes the artifacts of the previous phase and refuses to proceed if the upstream artifact is missing.

| # | Phase | Goal | Skill | Artifact produced |
|:--|:------|:-----|:------|:------------------|
| **1** | **Preflight** | Validate FASTQ + reference, decide `--platform`/`--size`/`--rounds`/`--taxon` and the execution environment | `bettamt-preflight/` | `params.json`, `params.rationale.md`, `nextflow.config` override |
| **2** | **Run** | Execute the workflow end-to-end with the override applied | `bettamt-run/` | Nextflow outputs (`results/**`) |
| **3** | **Debug** (on failure) | Interpret `.nextflow.log` / `work/` against the signature library | `bettamt-debug/` | `diagnosis.md` |
| **4** | **Assembly QC** | Completeness, circularity, NUMTs, publishability verdicts | `bettamt-qc/` | `report.md` |
| **5** | **Annotation QC** (optional, `--ref_gff`/`--ref_gb`) | Verify gene count, names, length tolerance, D-loop detection from `blast_extract_cds.py` output | `bettamt-annotate-qc/` | `annotation-report.md` |

> [!TIP]
> Read each skill's `SKILL.md` for the full procedure and its signature library (stderr pattern → cause → fix).

### Where the agent writes outputs

By default the agent operates in `$BETA_MT_HOME` (the pipeline directory). All handoff files are written *next to* the pipeline, not inside it, so the pipeline working directory stays clean. Set `RUN_DIR` to write them elsewhere — all five skills honor both variables.

### The handoff contract

| Direction | Artifact | Schema |
|---|---|---|
| Agent → Pipeline | `params.json` | Nextflow `-params-file` JSON; `params.rationale.md` explains every choice |
| Pipeline → Agent | `.nextflow.log`, `work/*/.command.*`, `results/**` | Nextflow-native files; the agent reads them, doesn't parse `.nf` internals |
| Agent → Human | `diagnosis.md`, `report.md`, `annotation-report.md` | Top-3-hypotheses / pass-warn-fail tables |

> The agent reasons; the pipeline executes.

## 🧰 Tools

All tools are resolved from conda-forge/bioconda via the pinned `BettaMt/pixi.toml`.

| Tool | Version | Role |
|:-----|:--------|:-----|
| `flye` | * | Long-read (ONT) mitogenome assembler |
| `getorganelle` | ≥1.7.7.1, <2 | Organelle/mitogenome assembly (Illumina path) |
| `circlator` | * | Assembly circularization |
| `racon` | * | Long-read polishing |
| `pypolca` | ≥0.4.0 | Illumina polishing (POLCA; pulls `bwa`, `freebayes`) |
| `bwa` | ≥0.7.17 | Alignment for the polca step |
| `minimap2` | * | ONT read alignment to reference |
| `samtools` | ≥1.18 | BAM/sam operations |
| `bbmap` | * | Read QC utilities |
| `seqkit` | * | Sequence stats and manipulation |
| `trnascan-se` | * | tRNA detection (`--taxon` model) |
| `bedtools` | * | Interval operations for annotation QC |
| `blast` | ≥2.14 | Gene annotation transfer (`blast_extract_cds.py`) |
| `python` | 3.11.6.* | Skill and helper scripts |

## 📁 Repository layout

```text
BettaMt-agents/
├── .agents/skills/            # the agentic layer (Agent Skills standard)
│   ├── bettamt-run/           # end-to-end orchestrator
│   ├── bettamt-preflight/     # params preparation
│   ├── bettamt-debug/         # failure interpreter
│   ├── bettamt-qc/            # post-assembly QC
│   └── bettamt-annotate-qc/   # post-annotation QC
├── BettaMt/                   # the pipeline (Nextflow DSL2, pixi)
│   ├── main.nf
│   ├── bin/
│   │   └── blast_extract_cds.py   # gene annotation helper
│   ├── conf/                  # base/local/slurm configs
│   ├── modules/
│   ├── subworkflows/          # ont_mitogenome.nf, illumina_mitogenome.nf
│   ├── pixi.toml              # pinned tool environment (pixi install)
│   └── README.md              # full parameter tables
├── test_data/
├── CITATION.cff               # rendered by GitHub as "Cite this repository"
└── README.md                  # this file
```

## 🔒 Hard guarantees

- **Filesystem evidence chain** — the boundary between agent and pipeline is the filesystem (`params.json` in, logs and `results/` out), not agent memory; the pipeline never reads the Markdown handoffs, the agent never parses `.nf` internals.
- **Reproducible environments** — every tool is pinned in `pixi.toml` and resolved via pixi; each Nextflow process invokes its tools through that manifest.
- **Rationale, not just output** — `params.rationale.md` records why every parameter was chosen, so a run's setup is auditable after the fact.
- **Bounded failure interpretation** — `bettamt-debug` maps observed stderr onto a signature library and writes ranked hypotheses; it does not free-restyle the log.
- **Clean working directory** — handoff files land next to the pipeline (`or `$RUN_DIR`), never inside it.

## Provenance

The five skills follow the [Agent Skills standard](https://agentskills.io/specification) and the AiX-BIO skill convention; this README follows the [bioinfo-skill-creator](https://github.com/cheahhl814/bioinfo-skill-creator) v1.2.0 README standard.

Upstream citations (Nextflow, Flye, GetOrganelle, Circlator, Racon, POLCA, tRNAscan-SE, `bait_mito`) are listed in [`CITATION.cff`](CITATION.cff), which GitHub renders as the "Cite this repository" button.

**License**: not yet set (`LicenseNotSet` in `CITATION.cff`) — default copyright applies; add a `LICENSE` (MIT / Apache-2.0 / GPL-3.0) before pointing external users at the repo.