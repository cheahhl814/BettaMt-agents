---
name: bettamt-qc
description: Phase 3 of the bettamt skill bundle — quality-check a finished BettaMt mitogenome (assembly QC: contiguity, circularity, tRNA set, coverage/NUMT risk, identity to reference; plus annotation QC for runs that used a reference GFF/GenBank). Works on outputs of BOTH runners — manual bash-recipe runs and the opt-in Nextflow runner. Use after a successful run, or when the user asks whether their mitogenome is complete, circular, NUMT-free, or publishable. Triggers: "bettamt QC", "check my mitogenome", "is the mitogenome complete", "post-assembly QC", "annotation QC", "bettamt report".
version: 1.0.0
updated: "2026-10-18"
triggers:
  - "bettamt QC"
  - "check my mitogenome"
  - "is the mitogenome complete"
  - "post-assembly QC"
  - "annotation QC"
  - "bettamt report"
---

# bettamt-qc — BettaMt post-assembly quality control

This skill runs **after** a successful run — manual bash-recipe or Nextflow mode. It never re-runs the runner. It produces one `report.md` with verdicts, evidence, and (optionally) a coverage profile plot; for annotated runs it also produces `annotation-report.md` (§6).

## 0. Locate the pipeline and the run

```bash
: "${BETA_MT_HOME:=$(cd "$SKILL_DIR/../.." && pwd)}"   # skill root = repo root
NF_HOME="$BETA_MT_HOME/runners/nextflow-runner"
test -f "$NF_HOME/main.nf" || { echo "Cannot locate the BettaMt runner tree. Set BETA_MT_HOME."; exit 1; }

# RUN_DIR is wherever the run landed (manual mode default: $BETA_MT_HOME/runs/<sample>)
RUN_DIR="${RUN_DIR:-$BETA_MT_HOME/runs/default}"
```

## 1. Identify platform and final contig

```bash
PARAMS="$RUN_DIR/params.json"
PLATFORM=$(jq -r '.platform' < "$PARAMS")

# Dual-mode output paths: manual recipes write run-root + stages/, NF writes results/
case "$PLATFORM" in
  ont)
    POLISHED=$(ls "$RUN_DIR"/*_mitogenome.fasta "$RUN_DIR"/results/polish/*_racon.fasta 2>/dev/null | head -1)
    TRIMMED_OR_BAITED="$(ls -d "$RUN_DIR"/stages "$RUN_DIR"/results/bait 2>/dev/null | head -1)"
    ;;
  illumina)
    POLISHED=$(ls "$RUN_DIR"/*_mitogenome.fasta "$RUN_DIR"/results/polish/*_polca.fasta 2>/dev/null | head -1)
    TRIMMED_OR_BAITED="$(ls -d "$RUN_DIR"/stages/report "$RUN_DIR"/results/qc 2>/dev/null | head -1)"
    ;;
  *) echo "Unknown platform: $PLATFORM"; exit 1 ;;
esac

test -n "$POLISHED" || { echo "Polished FASTA not found under results/polish/"; exit 1; }
```

## 2. Run the QC checks

Each check produces a verdict (✅ pass, ⚠️ warn, ❌ fail) and a one-line reason.

### Check 1 — contig size and contiguity

```bash
pixi run --manifest-path "$NF_HOME/pixi.toml" seqkit fx2tab -n -l "$POLISHED"
```

| Length                                   | Verdict                                               |
| ---------------------------------------- | ----------------------------------------------------- |
| 14 000–20 000 bp (typical vertebrate mt) | ✅                                                     |
| 10 000–14 000 or 20 000–25 000           | ⚠️ possible duplication or truncation                 |
| < 10 000 or > 25 000                     | ❌ likely not a mitogenome (or NUMT)                   |
| Multiple contigs of similar length       | ⚠️ fragmented — verify with `grep -c '>' "$POLISHED"` |

### Check 2 — circularity

```bash
# GetOrganelle outputs typically have a 75-bp terminal repeat (Hamming distance ~0)
pixi run --manifest-path "$NF_HOME/pixi.toml" seqkit head -n 2 "$POLISHED"  # first 2 records
# Compare first 75 bp of record 1 vs last 75 bp of last record
pixi run --manifest-path "$NF_HOME/pixi.toml" seqkit subseq -r 1:75 "$POLISHED"
# A circular mitogenome should show the same 75 bp at start and end
```

| Observation                                     | Verdict                                  |
| ----------------------------------------------- | ---------------------------------------- |
| First 75 bp == last 75 bp (exact or near-exact) | ✅ circularized                           |
| First 75 bp ≈ last 75 bp with 1–3 mismatches    | ✅ circular (Racon/POLCA polish artifact) |
| No obvious repeat                               | ⚠️ may be linear — check upstream contig |

### Check 3 — tRNA gene set

```bash
pixi run --manifest-path "$NF_HOME/pixi.toml" tRNAscan-SE -M vert -O -o /tmp/qc-trna.out "$POLISHED"
awk '$5 != "" {print $5}' /tmp/qc-trna.out | sort -u
```

| tRNAs found                           | Verdict                                                                  |
| ------------------------------------- | ------------------------------------------------------------------------ |
| 22 unique tRNAs (vertebrate standard) | ✅ complete                                                               |
| 20–21                                 | ⚠️ nearly complete; missing tRNA could be a real biology or a model miss |
| < 18                                  | ❌ incomplete — likely a fragment or wrong model                          |
| Run with `-M anonym` if < 15 found    | re-check; this should always pass after re-run                           |

### Check 4 — coverage profile (NUMT risk)

Map a downsampled subset of reads back to the polished contig and look for sudden coverage jumps (NUMT signature = a 5–10× local spike):

```bash
SAMPLE_READS=""
if [ "$PLATFORM" = "illumina" ]; then
    SAMPLE_READS=$(ls "$RUN_DIR"/results/qc/trimmed_1.fq.gz 2>/dev/null)
else
    SAMPLE_READS=$(ls "$RUN_DIR"/results/bait/*_mapped.fastq 2>/dev/null)
fi

pixi run --manifest-path "$NF_HOME/pixi.toml" minimap2 -ax $([ "$PLATFORM" = "illumina" ] && echo sr || echo map-ont) -t 4 "$POLISHED" "$SAMPLE_READS" \
  | pixi run --manifest-path "$NF_HOME/pixi.toml" samtools view -bF 4 - \
  | pixi run --manifest-path "$NF_HOME/pixi.toml" samtools sort -o /tmp/qc.bam -
pixi run --manifest-path "$NF_HOME/pixi.toml" samtools index /tmp/qc.bam
pixi run --manifest-path "$NF_HOME/pixi.toml" samtools depth -a /tmp/qc.bam > /tmp/qc.depth

# Summary stats
awk '{s+=$3; n++; if($3>max)max=$3; covs[n]=$3} END {
  print "mean=" s/n, "max=" max, "n_windows=" n
}' /tmp/qc.depth

# Flag a NUMT-suspicious window: median coverage × 5
awk '{a[NR]=$3} END {
  asort(a); med = (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2)
  if (max > 5*med) print "NUMT_SUSPICIOUS: max=" max " median=" med
}' /tmp/qc.depth
```

| Observation                 | Verdict                                       |
| --------------------------- | --------------------------------------------- |
| max coverage ≤ 5× median    | ✅ clean                                       |
| 5–10× median, single window | ⚠️ possible NUMT — inspect with IGV           |
| > 10× median, single window | ❌ strong NUMT signal — investigate that locus |

If the user wants a plot, generate a simple ASCII histogram of the depth file (or a PNG with `python -c "import matplotlib..."` if `matplotlib` is available — it is **not** in the BettaMt pixi env by default, so prefer ASCII).

### Check 5 — identity to reference

```bash
pixi run --manifest-path "$NF_HOME/pixi.toml" minimap2 -ax $([ "$PLATFORM" = "illumina" ] && echo sr || echo map-ont) "$REF_MITO" "$POLISHED" \
  | pixi run --manifest-path "$NF_HOME/pixi.toml" samtools view - \
  | awk '$6 != "*" { n++; m+=length($10)*($1/$1); } END { print "ref_pct_identity=" m/(n*length($10)) }'
# Simpler: just report alignment length and number of substitutions from CIGAR
```

| Identity to ref | Verdict                                     |
| --------------- | ------------------------------------------- |
| > 95%           | ✅ expected (same / closely related species) |
| 85–95%          | ⚠️ congeneric — verify expected divergence  |
| < 85%           | ❌ suspicious — might be a different locus   |

### Check 6 — GC content and skew

```bash
pixi run --manifest-path "$NF_HOME/pixi.toml" seqkit fx2tab -n -l -g "$POLISHED"
```

| GC%                       | Verdict                                                 |
| ------------------------- | ------------------------------------------------------- |
| 35–50% (most vertebrates) | ✅                                                       |
| 25–35 or 50–60%           | ⚠️ unusual but possible (e.g. AT-rich fish mitogenomes) |
| < 25% or > 60%            | ❌ suspicious                                            |

## 3. Output contract — write `report.md`

Always write to `$RUN_DIR/report.md` (and `$RUN_DIR/results/qc/` for archival). Template:

```markdown
# BettaMt assembly QC report

Sample:     $SAMPLE_ID
Platform:   $PLATFORM
Pipeline:   $BETA_MT_HOME
Contig:     $POLISHED
Generated:  <ISO8601>

## Verdict summary

| Check | Result | Notes |
|---|---|---|
| Contig size | ✅ / ⚠️ / ❌ | <length> bp; <one-line reason> |
| Circularity | ✅ / ⚠️ / ❌ | <evidence> |
| tRNA set    | ✅ / ⚠️ / ❌ | <N> tRNAs found (<list missing>) |
| Coverage profile | ✅ / ⚠️ / ❌ | mean=…× max=…× (NUMT risk <yes/no>) |
| Identity to ref  | ✅ / ⚠️ / ❌ | <pct>% over <N> bp |
| GC content       | ✅ / ⚠️ / ❌ | <pct>% |

**Overall**: <PASS / PASS-WITH-WARNINGS / FAIL>

## Evidence

### Contig size
\`\`\`
<raw seqkit output>
\`\`\`

### tRNA gene set
Found: <comma-separated list>
Missing: <comma-separated list, or "none">

### Coverage profile (ASCII histogram)
\`\`\`
<20-bin ASCII histogram of /tmp/qc.depth>
\`\`\`

### Identity to reference
<output of minimap2 alignment>

## Notes & recommended follow-ups
- <bullet list of any warnings or failures, with concrete next steps>

## Reproducibility
- params.json: $RUN_DIR/params.json
- pixi.lock:   $BETA_MT_HOME/pixi.lock
- this report: $RUN_DIR/report.md (regenerable via `bettamt-qc` skill)
```

If the user wants a PNG, generate one with:

```bash
pixi run --manifest-path "$NF_HOME/pixi.toml" python -c "
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
xs, ys = [], []
for line in open('/tmp/qc.depth'):
    p, d = line.split('\t')
    xs.append(int(p)); ys.append(int(d))
plt.figure(figsize=(10,3))
plt.plot(xs, ys)
plt.xlabel('position (bp)'); plt.ylabel('coverage')
plt.title('$SAMPLE_ID coverage profile')
plt.tight_layout(); plt.savefig('$RUN_DIR/coverage.png', dpi=150)
"
```

This requires adding `matplotlib` to `pixi.toml` if the user wants plots regularly — flag that as an opt-in.

## 4. What NOT to do

- Do **not** re-run the pipeline to "fix" issues — produce a report, let the user decide
- Do **not** fail the QC on a single check missing — flag it and explain
- Do **not** claim "NUMT" with high confidence unless coverage jump is > 10× median **and** the locus is supported by an alignment breakpoint
- Do **not** write to `$BETA_MT_HOME/...` paths — only `$RUN_DIR/...`

## 5. When to invoke other skills

- If QC fails because `--taxon` was wrong (Check 3 misses tRNAs) → recommend the user re-invoke `bettamt-preflight` with corrected params, then `nextflow run -resume`
- If coverage profile shows NUMT and it's a recurring problem → suggest sub-sampling the input reads at the `bait_mito` / `qc_short` stage, not at assembly
- If the report is "FAIL" on a non-recoverable axis (e.g. contig too short) → the user should re-design the experiment, not tweak params
- After this skill reports PASS and the run included annotation (`--ref_gff`/`--ref_gb`), run §6 below (annotation QC — folded from the former `bettamt-annotate-qc` skill). Annotation QC is meaningless on a broken assembly, so only run it after assembly QC has passed.

## 6. Annotation QC (folded from bettamt-annotate-qc)

Run only after the assembly QC verdict is ≥ PASS. Annotation outputs come from
`annotate_genes.sh` (manual) or the `ANNOTATE_GENES` process (nextflow mode) —
both wrap `bin/blast_extract_cds.py`. Locate them:

```bash
# manual mode: next to the polished fasta in the run stages; NF mode: results/annotation
ANNO_DIR="$(ls -d "$RUN_DIR"/annotated_genes.* 2>/dev/null && dirname "$(ls "$RUN_DIR"/annotated_genes.* | head -1)")"
[ -n "$ANNO_DIR" ] || ANNO_DIR="$RUN_DIR/results/annotation"
for f in annotated_genes.fasta annotated_genes.gff3 annotated_genes.bed; do
  test -f "$ANNO_DIR/$f" || echo "MISSING: $f"
done
```

### Check 1 — gene count

```bash
GENES=$(grep -c '^>' "$ANNO_DIR/annotated_genes.fasta")
echo "gene count: $GENES"
```

| Count | Verdict | Notes |
|---|---|---|
| 16 (13 PCG + 2 rRNA + D-loop) — `--no-trna` runs | ✅ | expected without tRNA annotation |
| 38 (16 + 22 tRNA) | ✅ | full vertebrate mitochondrial set |
| other | ⚠️ | list missing genes (compare `grep '^>' annotated_genes.fasta` against the reference GFF/GB gene names); a partial reference GFF legitimately yields fewer genes (see debug signature library) |

### Check 2 — per-gene length tolerance vs the reference

```bash
# bed vs reference-GFF CDS length; report genes whose transferred length deviates > 20%
awk -F'\t' '{print $4, $3-$2}' "$ANNO_DIR/annotated_genes.bed" | sort | head -20
grep -v '^#' "$REF_GFF" | awk -F'\t' '$3 ~ /gene|CDS|rRNA/ {split($9,a,"Name="); print a[2], $5-$4}' | sort | head -20
```
| Deviation | Verdict |
|---|---|
| ≤ 20 % | ✅ consistent BLAST transfer |
| 20–50 % | ⚠️ verify the gene manually (truncation/indel at the locus) |
| > 50 % | ❌ mis-annotation — check the debug signature library |

### Check 3 — D-loop detection

`blast_extract_cds.py` treats the largest inter-gene gap (≥ 100 bp) as the
control region / D-loop. If it reported `No D-loop detected`, check whether the
polished contig actually spans the control region (final length should be
14–20 kb; a truncated contig loses the D-loop first).

### Output — `annotation-report.md`

Write a per-gene verdict table (gene, length, deviation, verdict) plus counts,
and note the reference annotation used (GFF3 or GenBank, source species).
