---
name: bettamt-preflight
description: Prepare a Nextflow params.json for BettaMt from FASTQ + reference, plus a machine-specific nextflow.config override. Use when the user has raw reads and a reference mitogenome but has not yet decided --platform / --size / --rounds / --taxon / execution-environment. Triggers: "set up BettaMt", "preflight BettaMt", "what params for BettaMt", "BettaMt params.json", "bettamt preflight", "nextflow config local machine", "scale down BettaMt for laptop", "fix slurm partition BettaMt".
---

# bettamt-preflight — BettaMt parameter preparation

The pipeline `nextflow run main.nf` consumes a `-params-file params.json`. This skill produces that file plus a `params.rationale.md` audit trail. It does **not** run the pipeline.

## 0. Locate the pipeline

Resolve `BETA_MT_HOME`:
- Prefer `$BETA_MT_HOME` if set
- Otherwise look for a sibling directory named `BettaMt` next to this skill's parent (`../BettaMt`)
- Otherwise ask the user

```bash
: "${BETA_MT_HOME:=$(cd "$SKILL_DIR/../../../BettaMt" 2>/dev/null && pwd)}"
test -f "$BETA_MT_HOME/main.nf" || { echo "Cannot locate BettaMt. Set BETA_MT_HOME."; exit 1; }
```

All `pixi run` invocations below use `--manifest-path "$BETA_MT_HOME/pixi.toml"`.

## 1. Gather inputs from the user

Ask once, in one question batch if possible:

| Input | Required? | Default if absent |
|---|---|---|
| FASTQ path(s) — R1 for Illumina, single file for ONT | yes | — |
| R2 path (Illumina only) | conditional | — |
| `--ref_mito` reference mitogenome | yes | — |
| Reference annotation (GFF3 + FASTA, or GenBank) for gene annotation | no | skip `ANNOTATE_GENES` step |
| Organism / common name | no | `unknown` |
| Library strategy (WGS / mitogenome-capture / genome-skim) | no | `unknown` |
| Sequencer (ONT PromethION / MinION / Illumina NovaSeq / etc.) | no | infer from FASTQ |

## 2. Compute evidence (always run; never skip)

All commands use `pixi run` against the BettaMt env so versions are pinned.

### a) Read statistics

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" seqkit stats R1.fq.gz [R2.fq.gz]
```

Capture: `num_seqs`, `sum_len`, `avg_len`, `N50` (compute with `seqkit fx2tab -n -l R1.fq.gz | sort -k2 -n -r | awk ...` if not printed).

### b) Reference length and GC

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" seqkit fx2tab -n -l -g ref_mito.fasta
```

### c) Coverage estimate

For Illumina paired-end (use first 100 k read pairs to keep it fast):

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" minimap2 -ax sr -t 4 ref_mito.fasta <(seqkit head -n 100000 R1.fq.gz) <(seqkit head -n 100000 R2.fq.gz) \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools view -bF 4 - \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools depth -a - \
  | awk '{s+=$3; n++} END {print "mean_cov=" s/n}'
```

For ONT:

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" minimap2 -ax map-ont -t 4 ref_mito.fasta <(seqkit head -n 50000 R1.fq.gz) \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools view -bF 4 - \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools depth -a - \
  | awk '{s+=$3; n++} END {print "mean_cov=" s/n}'
```

### d) tRNA content of the reference (drives `--taxon`)

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" tRNAscan-SE -M vert -O -o /tmp/trna.out ref_mito.fasta
awk '$5 == "Phe" || $5 == "His" || $5 == "Pro"' /tmp/trna.out
```

If Phe/His/Pro are all present → `vertebrate` model is appropriate. If absent → fall back to `anonym` and warn.

### e) Platform auto-detection (if not told)

- Read lengths mean > 1000 bp → ONT
- Read lengths mean < 1000 bp AND two FASTQ files → Illumina
- If ambiguous, ask.

### f) Disk space check

BettaMt is small compared to nf-betta (mitogenome assembly, not whole genome), but it can still touch 50–200 GB on a high-coverage Illumina run (intermediates from SPAdes). Run the same pre-check as the pipeline will need:

```bash
df -BG "$RUN_DIR" | awk 'NR==2 {print "Free space:", $4, "(recommend > 50 GB for typical mitogenome assembly, > 200 GB for high-coverage Illumina)"}'
```

If free space < 50 GB, warn — intermediates are mostly in `work/` and can be deleted after the run finishes, but during the run they coexist with `results/`.

## 2.5 Profile and machine-spec detection (run before §3)

BettaMt ships with three Nextflow config files (`nextflow.config` + `conf/base.config` + `conf/{local,slurm}.config`):

- `local` — `executor = 'local'` only; per-process resources fall back to `base.config` defaults (`cpus = 1, memory = '4 GB'`). Intentionally minimal: a laptop can run BettaMt's mitogenome assembly in tens of minutes to a few hours without explicit tuning.
- `slurm` — per-process `withName:` blocks with hardcoded `clusterOptions` for partitions `cpu1.medm`, `cpu2.largem`, `cpu1`, `cpu2`. The two heaviest processes are `BAIT_MITO` (24 cpus / 180 GB) and `GETORGANELLE_ASM` (16 cpus / 96 GB).
- `base.config` — applied to every executor; sets `cpus = 1`, `memory = '4 GB'`, and the retry policy.

The base-config defaults mean **a `local` user gets 1 cpu / 4 GB per process** unless they explicitly override. For a mitogenome that's often enough, but the heavier processes (GetOrganelle, Flye) can blow past 4 GB. This section auto-detects the environment, compares the machine to what the pipeline expects, and writes a user-side override `nextflow.config.<env>.local` that you apply at run time via `-c`. The original config files are never edited.

### a) Ask the execution environment (one question)

```
? Execution environment
  - slurm    : HPC cluster with SLURM scheduler
  - local    : workstation (no scheduler)
  - docker   : Docker / Singularity container (auto-detected if $SINGULARITY_CONTAINER or $CONTAINER env is set)
```

If the user picks `slurm`, also ask which partitions are valid on their cluster (or run `sinfo` to list them; see §2.5d). If they pick `local`, run the resource check in §2.5c to decide whether to add per-process overrides.

### b) Detect environment (auto; never skip)

```bash
# Are we already on a SLURM cluster?
if [ -n "$SLURM_JOB_ID" ] || [ -n "$SLURM_CLUSTER_NAME" ]; then
  ENV_HINT=slurm
elif [ -n "$SINGULARITY_CONTAINER" ] || [ -n "$APPTAINER_CONTAINER" ] || [ -n "$CONTAINER" ]; then
  ENV_HINT=container
else
  ENV_HINT=local
fi

echo "Detected execution environment: $ENV_HINT"
```

This is a hint, not a final choice — the user’s answer in §2.5a wins. If their answer disagrees with the hint, surface it (“you picked slurm, but no $SLURM_* env vars are set — are you submitting from a login node?”) and proceed with the user’s choice.

### c) Resource check (for `local` environment)

```bash
NCPU=$(nproc)
RAM_GB=$(free -g | awk '/^Mem:/{print $2}')
DISK_GB=$(df -BG "$RUN_DIR" | awk 'NR==2 {gsub("G",""); print $2}')
echo "Local machine: ncpu=$NCPU, ram=${RAM_GB}G, disk=${DISK_GB}G free"
echo ""

# Read the slurm profile's per-process maxima as the upper bound. We use
# slurm's values because the local profile has none (defaults to base.config's
# 1 cpu / 4 GB), so anything you'd run on slurm is also what you'd want
# locally if the machine can handle it.
SLURM_LIMITS=$(awk '/^    slurm \{/,/^    \}/' "$BETA_MT_HOME/conf/slurm.config" \
  | awk '/withName:/{name=$2; gsub(/['\''{}]/,"",name)} /cpus   =/{cpu[name]=$3} /memory =/{gsub(/'\''/,"",$3); sub(/GB/,"",$3); mem[name]=$3} END {for (k in cpu) printf "%s\t%d\t%d\n", k, cpu[k], mem[k]}')
echo "Pipeline's per-process maxima (from conf/slurm.config):"
echo "$SLURM_LIMITS" | column -t -s $'\t'
echo ""
echo "Note: conf/base.config defaults to cpus=1, memory=4 GB for all"
echo "processes. The local profile inherits this. Without an override,"
echo "each process gets exactly that — which is fine for most mitogenome"
echo "jobs but will OOM on GetOrganelle (SPAdes) and Flye."
```

**Scaling rules for `local`** (apply per process name):

| Process | Default (base.config) | Local override when host has ≥ 16 GB RAM |
| --- | --- | --- |
| `BAIT_MITO` (ONT)         | 1 cpu / 4 GB | `cpus = min(24, ncpu)`, `memory = '16 GB'` (bait mapping is light) |
| `FLYE` (ONT)              | 1 cpu / 4 GB | `cpus = min(12, ncpu)`, `memory = min(64, floor(ram_gb * 0.6)) GB` |
| `GETORGANELLE_FILTER`     | 1 cpu / 4 GB | `cpus = min(6, ncpu)`, `memory = '8 GB'` |
| `GETORGANELLE_ASM`        | 1 cpu / 4 GB | `cpus = min(16, ncpu)`, `memory = min(96, floor(ram_gb * 0.7)) GB` (SPAdes is the heaviest) |
| `RACON` (ONT polish)      | 1 cpu / 4 GB | `cpus = min(6, ncpu)`, `memory = min(24, floor(ram_gb * 0.3)) GB` |
| `POLCA` (Illumina polish) | 1 cpu / 4 GB | `cpus = min(8, ncpu)`, `memory = min(32, floor(ram_gb * 0.4)) GB` |
| `CIRCLATOR`               | 1 cpu / 4 GB | `cpus = 4, memory = '8 GB'` |
| `FIRSTGENE`               | 1 cpu / 4 GB | `cpus = 4, memory = '8 GB'` |
| `QC_SHORT`                | 1 cpu / 4 GB | `cpus = 4, memory = '8 GB'` |

**`memory` is hard**: if a process needs `X GB` and the node has only `Y < X`, the process gets OOM-killed regardless of `cpus`. So if the user’s RAM is short, the override should *reduce memory requests* and accept longer walltime, or run on SLURM.

If the user’s machine can’t meet the heavier processes (`GETORGANELLE_ASM` ≥ 96 GB, `BAIT_MITO` ≥ 180 GB), the skill should:

1. **Warn loudly** — these processes may fail or swap-thrash.
2. **Suggest a fallback** — either accept the reduced-memory override and longer walltime, or push the user to SLURM.
3. **Offer to generate the override anyway** — the user may have swap or a beefy hidden node.

### d) SLURM partition check (for `slurm` environment)

```bash
# What partitions are actually available?
if command -v sinfo >/dev/null 2>&1; then
  AVAILABLE_PARTITIONS=$(sinfo -h -o "%P" | sort -u | grep -v "^$" | tr '\n' ' ')
  echo "Partitions on this cluster: $AVAILABLE_PARTITIONS"
else
  echo "sinfo not available — cannot check partitions; will fall back to user-supplied names."
  AVAILABLE_PARTITIONS=""
fi

# What does conf/slurm.config hardcode?
HARDCODED_PARTITIONS=$(grep -oE "partition=[^ ]+" "$BETA_MT_HOME/conf/slurm.config" | sort -u)
echo "Hardcoded in conf/slurm.config: $HARDCODED_PARTITIONS"

# Find hardcoded partitions that don't exist on this cluster
if [ -n "$AVAILABLE_PARTITIONS" ]; then
  for hp in $HARDCODED_PARTITIONS; do
    name=${hp#partition=}
    if ! echo " $AVAILABLE_PARTITIONS " | grep -q " $name "; then
      echo "  WARNING: partition '$name' hardcoded in conf/slurm.config is NOT available here"
    fi
  done
fi
```

If any hardcoded partition is missing, prompt the user for the correct partition name on their cluster and write it into the override. Note that BettaMt's `clusterOptions` also include `--job-name` and `-o/-e` paths; you typically only need to override the partition flag.

### e) Write the user-side override file

Generate `nextflow.config.<env>.local` (where `<env>` = `local` or `slurm` or `container`) **next to `params.json`** in the run directory. The user applies it via `-c`:

```bash
# Example for local environment on a 16-core / 64 GB workstation:
cat > nextflow.config.local.local <<'EOF'
// Auto-generated by bettamt-preflight on <ISO8601>
// Detected: nproc=16, ram=64G, disk=500G
// Applied at run time:  nextflow run main.nf -profile local -c nextflow.config.local.local -params-file params.json

profiles {
    local {
        process {
            withName: 'BAIT_MITO'           { cpus = 16 ; memory = '40 GB'  }
            withName: 'FLYE'                { cpus = 12 ; memory = '40 GB'  }
            withName: 'GETORGANELLE_ASM'    { cpus = 8  ; memory = '48 GB'  }
            withName: 'GETORGANELLE_FILTER' { cpus = 6  ; memory = '8 GB'   }
            withName: 'RACON'               { cpus = 6  ; memory = '16 GB'  }
            withName: 'POLCA'               { cpus = 8  ; memory = '24 GB'  }
            withName: 'CIRCLATOR'           { cpus = 4  ; memory = '8 GB'   }
            withName: 'FIRSTGENE'           { cpus = 4  ; memory = '8 GB'   }
            withName: 'QC_SHORT'            { cpus = 4  ; memory = '8 GB'   }
        }
    }
}
EOF
```

**Why `withName:` and quoted process names?** BettaMt targets processes by their DSL2 process name (`'BAIT_MITO'`, `'FLYE'`, etc.), not by generic labels. The quotes are required because process names are strings in Nextflow config syntax. Do **not** use `withLabel:` here — it will silently do nothing because BettaMt does not set `label` directives on its processes.

The original `conf/slurm.config` and `conf/base.config` are **never modified**. The override is loaded via `-c nextflow.config.local.local` at run time. If the user re-runs on a different machine, they can re-generate the override or edit it.

**Naming convention** (the file lives next to `params.json`, not inside `$BETA_MT_HOME`):
- `nextflow.config.local.local`   — workstation override
- `nextflow.config.slurm.local`   — SLURM cluster override
- `nextflow.config.container.local` — container override

The `nextflow.config.<env>.local` name is intentional: `<env>` matches the `-profile <env>` so the user can run `nextflow run main.nf -profile <env> -c nextflow.config.<env>.local ...` without thinking.

### f) What if the user wants to skip the override?

If the user’s machine has plenty of RAM (≥ 180 GB for ONT or ≥ 96 GB for Illumina + ≥ 16 cores), the slurm profile’s per-process values can be used directly. For `local` they still need an override because the `local` profile inherits only the 1-cpu / 4-GB base defaults. The skill should *still write the override file* as a fallback — it costs nothing and is reusable. But the run-command template in §4 should show the `-c` flag as the recommended form.

If the user explicitly says “I don’t want an override”, honor that. The skill then just emits `params.json` + `params.rationale.md` and leaves the user to pick `-profile local` or `-profile slurm` raw.

## 3. Reasoning framework

Apply these rules. Each recommendation **must** cite the evidence in `params.rationale.md`.

### `--platform`

| Evidence | Recommendation |
|---|---|
| User says "Illumina" / "NovaSeq" / "paired" | `illumina` |
| User says "ONT" / "PromethION" / "MinION" / "long-read" | `ont` |
| mean read len > 1000 bp (computed) | `ont` |
| two FASTQ files provided | `illumina` |
| All else fails | **ask the user** |

### `--size` (ONT only)

- Default: `ref_len × 1.05` rounded to nearest `k` (e.g. `16k` for 16 200 bp)
- If `ref_len < 5 000`: warn "reference looks like a fragment, expect fragmented assembly"
- If `ref_len > 25 000`: warn "reference is longer than typical mitogenome — likely nuclear scaffold or NUMT"

### `--rounds` (Illumina only)

| Coverage (Illumina) | `--rounds` |
|---|---|
| < 50× | 30 |
| 50–500× | 15 (default) |
| 500–5 000× | 20 |
| > 5 000× (very high, e.g. WGS) | 30 + warn "GetOrganelle can produce tangled graphs at very high coverage; consider `--reduce-reads-for-coverage 2000` via the module override" |

### `--taxon`

- tRNAscan finds Phe/His/Pro with the reference → `vertebrate`
- tRNAscan finds nothing meaningful → `anonym` + warn "could not auto-detect taxon; consider providing a closer reference"
- User explicitly says "mammal" / "invertebrate" → respect it

### `--sample_id` (Illumina only)

- Default: strip `_R1`/`_R2`/`.fastq`/`.fq`/`.gz` from the R1 basename
- Example: `betta_R1.fastq.gz` → `betta`

### Gene annotation parameters (optional, both platforms)

If the user supplied a reference annotation, propagate it into `params.json`. Accept **either** GFF3 + companion FASTA **or** a single GenBank file — never both.

- `--ref_gff` and `--ref_gff_fasta`: GFF3 + companion reference genome FASTA
- `--ref_gb`: self-contained GenBank file

If the user did **not** supply a reference annotation, set all three keys to `""` in `params.json` and `ANNOTATE_GENES` will be skipped silently. Never invent one.

If the user is unsure whether to provide one, default to skipping. Annotation is most useful when the reference is the same species or a close congener; a divergent reference produces noisy BLAST hits that fail the QC length-tolerance check.

## 4. Output contract — write exactly these files

### `params.json` (hard schema, all keys required unless marked optional)

```json
{
  "platform":  "ont|illumina",
  "reads":     "/abs/path/to/R1_or_ONT.fastq.gz",
  "reads_r2":  "/abs/path/to/R2.fastq.gz",
  "ref_mito":  "/abs/path/to/reference.fasta",
  "size":      "16k",
  "taxon":     "vertebrate",
  "rounds":    15,
  "sample_id": "my_sample",
  "ref_gff":       "",
  "ref_gff_fasta": "",
  "ref_gb":        "",
  "trnascan":      "",
  "skip_trna":     false
}
```

- For ONT: omit `reads_r2` and `sample_id` (or set to `null` and the dispatcher will ignore)
- For Illumina: `rounds` is required; `size` can be omitted
- `taxon` is always required (shared by both paths)
- All five annotation keys (`ref_gff`, `ref_gff_fasta`, `ref_gb`, `trnascan`, `skip_trna`) are required; leave `""`/`false` to skip annotation

### `nextflow.config.<env>.local` (machine-specific override — see §2.5)

This file lives next to `params.json` in the run directory and is applied via `-c`. Filename format: `nextflow.config.<env>.local` where `<env>` matches the chosen `-profile`. **Never edit `$BETA_MT_HOME/conf/{base,local,slurm}.config` in place** — always use a user-side override.

See §2.5e for the template. The contents depend on the detected machine; do not write a fixed stub.

### `params.rationale.md` (audit trail — **do not skip**)

A short markdown file with one bullet per decision. Template:

```markdown
# BettaMt preflight rationale

Generated: <ISO8601>
Pipeline:  $BETA_MT_HOME (commit <hash if available>)

## Evidence
- reads:     R1=…R1.fq.gz, R2=…R2.fq.gz, n=… × 2, mean_len=… bp, N50=… bp
- reference: ref.fasta, len=… bp, GC=…%
- coverage:  …× (from downsampled minimap2)
- tRNAscan:  Phe=yes, His=yes, Pro=yes (vertebrate marker set present)
- disk:      …G free (recommend > 50 GB; > 200 GB for high-coverage Illumina)

## Environment
- env_hint (auto):   <slurm|local|container>
- env (user):        <slurm|local|container>
- machine:           nproc=…, ram=…G, disk=…G free
- slurm_partitions:  <list from sinfo, or "not on slurm">
- hardcoded_partitions_in_config:  <list from conf/slurm.config>
- partition_mismatch:               <none|<partitions not on this cluster>>

## Config override
- wrote:        nextflow.config.<env>.local
- scaling_rule: <"none (machine matches profile)" | "scaled per-process via withName: blocks; cpus=min(config, nproc); memory=floor(ram_gb * fraction)">
- warnings:     <list, e.g. "GETORGANELLE_ASM memory 96 GB > local ram 64 GB — expect OOM; suggest reducing coverage or running on slurm">

## Decisions
- --platform ont / illumina  : <reason>
- --size 16k (ONT)          : ref_len=… × 1.05 = … → rounded to 16k
- --rounds 15 (Illumina)    : coverage …× falls in 50–500× band
- --taxon vertebrate        : tRNAscan found Phe/His/Pro marker set
- --sample_id foo           : derived from R1 basename
- --ref_gff / --ref_gb      : <provided; same species / congener / skipped>

## Warnings
- (none, or list any that fired)

## Run command
\`\`\`bash
cd "$BETA_MT_HOME"
nextflow run main.nf -profile <env> -c nextflow.config.<env>.local -params-file params.json
# If you skipped the override (machine exactly matches the profile):
nextflow run main.nf -profile <env> -params-file params.json
\`\`\`
```

## 5. What NOT to do

- Do **not** invent parameters the schema doesn't define
- Do **not** guess `--taxon` without running `tRNAscan-SE` on the reference
- Do **not** run `nextflow run` — this skill produces the params, the user (or a follow-up skill) runs the pipeline
- Do **not** write to `$BETA_MT_HOME/results/` — that's the pipeline's own output directory
- Do **not** edit `$BETA_MT_HOME/conf/{base,local,slurm}.config` in place — always use a `nextflow.config.<env>.local` user override via `-c` (see §2.5)
- Do **not** use `withLabel:` in the BettaMt override — BettaMt targets processes by `withName:` (e.g. `'BAIT_MITO'`, `'FLYE'`, `'GETORGANELLE_ASM'`). A `withLabel:` block will silently do nothing.
- Do **not** recommend `-profile local` for high-coverage Illumina without checking that the host has ≥ 96 GB RAM — `GETORGANELLE_ASM` (SPAdes) will OOM on a 16-GB laptop
- Do **not** skip the partition-mismatch check on slurm — a hardcoded `--partition=cpu1.medm` that doesn't exist on the user's cluster will fail every task submission with a cryptic SLURM error
- Do **not** expect GPU support — BettaMt has no GPU-friendly process. If the user wants faster polishing, suggest `-profile slurm` with more CPU allocation, not GPU.

## 6. Handoff

After writing `params.json` + `params.rationale.md`, print the run command and stop. If the user says "go", invoke the run command yourself (with `--resume` for re-runs).
