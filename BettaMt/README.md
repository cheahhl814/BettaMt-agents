# BettaMt — Mitochondrial Genome Assembly Pipeline

A [Nextflow DSL2](https://www.nextflow.io/) pipeline for *de novo* assembly of mitochondrial genomes from either **Oxford Nanopore long reads** or **Illumina paired-end short reads**. It covers read QC, baiting, assembly, organelle extraction, circularization, rotation, and platform-appropriate polishing — all tools managed through a reproducible [pixi](https://pixi.sh) environment.

## Pipeline Overview

### Long-read path (`--platform ont`)

```
reads + ref_mito
      │
      ▼
 1. BAIT_MITO         — Filter reads by length and map to a reference mitogenome
      │
      ▼
 2. FLYE              — Assemble baited reads into contigs
      │
      ▼
 3. GETORGANELLE_FILTER — Extract the mitochondrial contig from the assembly graph
      │                       (also runs `circlator clean` internally)
      ├──────────────────────────────┐
      ▼                              ▼
 4. FIRSTGENE          5. CIRCLATOR fixstart — Rotate to trnF start site
    (trnF / Phe)                       │
      │                                ▼
      └─────────────────► 6. RACON  — Polish with ONT raw reads
```

### Short-read path (`--platform illumina`)

```
R1 + R2  (paired-end FASTQ)
      │
      ▼
 1. QC_SHORT          — sequali QC + fastp adapter/quality trim
      │
      ▼
 2. GETORGANELLE_ASM  — De-Bruijn assembly directly from reads
      │                  (uses --ref_mito as initial seed)
      ├──────────────────────────────┐
      ▼                              ▼
 3. FIRSTGENE         4. CIRCLATOR fixstart — Rotate to trnF start site
    (trnF / Phe)                       │
      │                                ▼
      └────────────────► 5. POLCA    — Polish with paired Illumina reads
                                         (--careful mode)
```

## Dependencies

All bioinformatics tools are managed by `pixi` and declared in `pixi.toml`. Only Nextflow itself needs to be installed separately.

| Tool | Role | Path |
|---|---|---|
| Nextflow | Workflow orchestration | external |
| Python ≥ 3.11 | Runtime for GetOrganelle / pypolca / tRNAscan-SE wrapper | both |
| SeqKit | ONT read length filtering | ONT |
| BBTools (bbduk) | k-mer–based read baiting | both |
| Minimap2 | Long-read alignment (bait + Racon polish) | ONT |
| Samtools | SAM/BAM processing (≥ 1.18 for POLCA) | both |
| Flye | Long-read de novo assembly | ONT |
| GetOrganelle | Organelle contig extraction (Flye post-filter) | ONT |
| GetOrganelle | Short-read organelle assembler | Illumina |
| Circlator | Contig cleaning and rotation | both |
| tRNAscan-SE | tRNA gene annotation | both |
| BEDTools | Sequence extraction from BED coordinates | both |
| Racon | Long-read polishing | ONT |
| **Sequali** | Per-base / per-read QC | Illumina |
| **fastp** | Adapter + quality trimming | Illumina |
| **BWA** | Short-read alignment (for POLCA) | Illumina |
| **freebayes** | Variant calling (for POLCA, via pypolca) | Illumina |
| **pypolca** | Short-read polishing (POLCA re-impl.) | Illumina |

## Installation

**1. Install Nextflow** (requires Java ≥ 11):
```bash
pixi global install -c conda-forge -c bioconda nextflow
```

**2. Install all pipeline tools** via the bundled pixi environment:
```bash
cd BettaMt/
pixi install
```

This resolves all tool versions declared in `pixi.toml` and locks them in `pixi.lock` for full reproducibility.

## Usage

### Long-read (ONT)

```bash
nextflow run main.nf -profile local \
  --platform ont \
  --reads    data/reads.fastq \
  --ref_mito data/reference.fasta \
  --size     16k \
  --taxon    vertebrate
```

| Parameter | Description | Example |
|---|---|---|
| `--reads` | Single ONT FASTQ (gzipped OK) | `data/reads.fastq.gz` |
| `--ref_mito` | Reference mitogenome (used for k-mer baiting) | `data/reference.fasta` |
| `--size` | Estimated mitogenome size | `16k` |
| `--taxon` | tRNAscan-SE model (`mammal` / `vertebrate`) | `vertebrate` |

### Short-read (Illumina)

```bash
nextflow run main.nf -profile local \
  --platform   illumina \
  --reads      data/sample_R1.fastq.gz \
  --reads_r2   data/sample_R2.fastq.gz \
  --ref_mito   data/reference.fasta \
  --rounds     15 \
  --sample_id  betta_splendens \
  --taxon      vertebrate
```

| Parameter | Description | Example |
|---|---|---|
| `--reads` | R1 paired-end FASTQ (gzipped OK) | `data/sample_R1.fastq.gz` |
| `--reads_r2` | R2 paired-end FASTQ (gzipped OK) | `data/sample_R2.fastq.gz` |
| `--ref_mito` | Reference mitogenome (used as GetOrganelle seed) | `data/reference.fasta` |
| `--rounds` | GetOrganelle `-R` (extension rounds; 10–15 typical) | `15` |
| `--sample_id` | Sample identifier (defaults to `--reads` basename) | `betta_splendens` |
| `--taxon` | tRNAscan-SE model (`mammal` / `vertebrate`) | `vertebrate` |

### Cluster execution

```bash
nextflow run main.nf -profile slurm \
  --platform   illumina \
  --reads      data/sample_R1.fastq.gz \
  --reads_r2   data/sample_R2.fastq.gz \
  --ref_mito   data/reference.fasta
```

## Output

Results are written to `results/` within the project directory:

```
results/
├── bait/            # ONT: length-filtered and mapped reads
├── qc/              # Illumina: sequali + fastp reports
├── <sample_id>/     # Illumina: GetOrganelle assembly workdir (one per sample)
├── flye/            # ONT: Flye assembly output (contigs, assembly graph)
├── circlator/       # Rotated, circularized contig (both platforms)
└── polish/          # *_racon.fasta (ONT) or *_polca.fasta (Illumina)
```

## Project layout

```
BettaMt/
├── main.nf                       # entry point, dispatches by --platform
├── nextflow.config               # profile plumbing
├── pixi.toml                     # tool versions
├── conf/
│   ├── base.config               # shared process defaults
│   ├── local.config              # local-executor profile
│   └── slurm.config              # SLURM per-process resources
├── modules/local/                # one process per file, UPPER_CASE
│   ├── bait_mito.nf
│   ├── flye.nf
│   ├── getorganelle_filter.nf
│   ├── firstgene.nf
│   ├── circlator.nf
│   ├── racon.nf
│   ├── qc_short.nf               # Illumina: sequali + fastp
│   ├── getorganelle_asm.nf       # Illumina: short-read assembler
│   └── polca.nf                  # Illumina: pypolca short-read polisher
└── subworkflows/
    ├── ont_mitogenome.nf
    └── illumina_mitogenome.nf
```

## Citations

* Kolmogorov M. *et al.* — **Flye** (Nat Biotechnol 2019)
* Jin J.-J. *et al.* — **GetOrganelle** (Mol Biol Evol 2020 / Genome Biol 2020)
* Hunt M. *et al.* — **Circlator** (Genome Biol 2015)
* Vaser R. *et al.* — **Racon** (Genome Res 2017)
* Chen S. *et al.* — **fastp** (Bioinformatics 2018)
* Zimin A.V. & Salzberg S.L. — **POLCA** (PLoS Comput Biol 2020)
* Bouras G. *et al.* — **pypolca** (Python re-implementation of POLCA)
* Li H. — **minimap2** / **BWA** (Bioinformatics 2016 / 2009)
