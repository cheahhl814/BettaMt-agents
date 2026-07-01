---
name: bettamt-preflight
description: Prepare a Nextflow params.json for BettaMt from FASTQ + reference. Use when the user has raw reads and a reference mitogenome but not yet decided --platform / --size / --rounds / --taxon. Triggers: "set up BettaMt", "preflight BettaMt", "what params for BettaMt", "BettaMt params.json", "bettamt preflight".
---

# bettamt-preflight — BettaMt parameter preparation

The pipeline `nextflow run main.nf` consumes a `-params-file params.json`. This skill produces that file plus a `params.rationale.md` audit trail. It does **not** run the pipeline.

## 0. Locate the pipeline

Resolve `BETA_MT_HOME`:
- Prefer `$BETA_MT_HOME` if set
- Otherwise look for a sibling directory named `BettaMt` next to this skill's parent (`../BettaMt`)
- Otherwise ask the user

```bash
: "${BETA_MT_HOME:=$(cd "$(dirname "$SKILL_DIR")/../BettaMt" 2>/dev/null && pwd)}"
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
  "sample_id": "my_sample"
}
```

- For ONT: omit `reads_r2` and `sample_id` (or set to `null` and the dispatcher will ignore)
- For Illumina: `rounds` is required; `size` can be omitted
- `taxon` is always required (shared by both paths)

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

## Decisions
- --platform ont / illumina  : <reason>
- --size 16k (ONT)          : ref_len=… × 1.05 = … → rounded to 16k
- --rounds 15 (Illumina)    : coverage …× falls in 50–500× band
- --taxon vertebrate        : tRNAscan found Phe/His/Pro marker set
- --sample_id foo           : derived from R1 basename

## Warnings
- (none, or list any that fired)

## Run command
\`\`\`bash
cd "$BETA_MT_HOME"
nextflow run main.nf -profile slurm -params-file params.json
\`\`\`
```

## 5. What NOT to do

- Do **not** invent parameters the schema doesn't define
- Do **not** guess `--taxon` without running `tRNAscan-SE` on the reference
- Do **not** run `nextflow run` — this skill produces the params, the user (or a follow-up skill) runs the pipeline
- Do **not** write to `$BETA_MT_HOME/results/` — that's the pipeline's own output directory

## 6. Handoff

After writing `params.json` + `params.rationale.md`, print the run command and stop. If the user says "go", invoke the run command yourself (with `--resume` for re-runs).
