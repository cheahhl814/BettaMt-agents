---
name: bettamt-annotate-qc
description: Quality-check the gene annotation output from blast_extract_cds.py. Use after a successful BettaMt run that included --ref_gff/--ref_gb. Verifies gene count, gene names, length-tolerance vs reference, and D-loop detection. Triggers: "annotate QC", "check annotation", "did the BLAST annotation work", "bettamt annotate report".
---

# bettamt-annotate-qc — post-annotation quality control

This skill runs **after** `bettamt-qc` has confirmed the assembly is sane. It never re-runs the pipeline. It produces one `annotation-report.md` with per-gene verdicts.

## 0. Locate the pipeline and the run

```bash
: "${BETA_MT_HOME:=$(cd "$(dirname "$SKILL_DIR")/../BettaMt" 2>/dev/null && pwd)}"
test -f "$BETA_MT_HOME/main.nf" || { echo "Cannot locate BettaMt. Set BETA_MT_HOME."; exit 1; }

RUN_DIR="${RUN_DIR:-$BETA_MT_HOME}"
ANNO_DIR="$RUN_DIR/results/annotation"
test -d "$ANNO_DIR" || { echo "No annotation output at $ANNO_DIR"; exit 1; }
```

## 1. Confirm annotation files exist

```bash
test -f "$ANNO_DIR/annotated_genes.fasta" || { echo "Missing annotated_genes.fasta"; exit 1; }
test -f "$ANNO_DIR/annotated_genes.gff3"  || { echo "Missing annotated_genes.gff3"; exit 1; }
test -f "$ANNO_DIR/annotated_genes.bed"   || { echo "Missing annotated_genes.bed"; exit 1; }
```

## 2. Run the annotation QC checks

### Check 1 — gene count

```bash
grep -c '^>' "$ANNO_DIR/annotated_genes.fasta"
```

| Count | Verdict | Notes |
|---|---|---|
| 16 (13 PCG + 2 rRNA + D-loop) when `--no-trna` | ✅ | expected for `--no-trna` runs |
| 38 (16 + 22 tRNA) | ✅ | expected full vertebrate mitochondrial set |
| 15 (missing D-loop) | ⚠️ | D-loop not detected — verify largest gap < 100 bp or `--target` not circularized |
| 16–37 (partial tRNAs) | ⚠️ | tRNAscan-SE didn't find all 22 tRNAs; possibly wrong `--taxon` model or truncated mitogenome |
| < 13 | ❌ | too few PCGs annotated — possible reference is too divergent or assembly is broken |
| 0 | ❌ | BLAST failed for all genes; check reference compatibility |

### Check 2 — gene name normalization

Extract gene names from the FASTA headers and confirm none are unnormalized aliases:

```bash
grep '^>' "$ANNO_DIR/annotated_genes.fasta" | sed 's/^>//' | sort -u
```

| Observation | Verdict |
|---|---|
| All names match: `ND1, ND2, COXI, COXII, COXIII, ATPase6, ATPase8, ND3, ND4, ND4L, ND5, ND6, Cytb, 12S, 16S, D-loop, tRNA-XXX` | ✅ |
| Names like `COX1`, `ATP6`, `CYTB`, `RNR1`, `RNR2` (un-aliased) appear | ❌ | the `GENE_ALIASES` dict in `blast_extract_cds.py` should have caught these; report a bug |
| Names like `gene-XXX`, `cds-XXX`, `YP_003024028.1` appear | ❌ | the parser used a protein ID instead of a gene name; report a bug |

### Check 3 — gene order and genomic span

```bash
cut -f1,2,3,4,6 "$ANNO_DIR/annotated_genes.bed" | head -30
```

Compute the total span of the annotation (sum of gene lengths):

```bash
awk '{s+=$3-$2} END {print "annotated_bp=" s}' "$ANNO_DIR/annotated_genes.bed"
```

| Observation | Verdict |
|---|---|
| annotated_bp ≈ 15,000–16,000 bp (13 PCGs + 2 rRNAs) | ✅ |
| annotated_bp ≈ 11,000–13,000 bp (missing some genes) | ⚠️ partial annotation |
| annotated_bp > 17,000 bp | ❌ over-annotation (e.g. duplicated genes or a spurious feature) |

The total annotated span should be roughly the polished contig length minus the D-loop. The contig length itself is in the `bettamt-qc` report; cross-check that this skill's annotated span + D-loop length ≈ contig length.

### Check 4 — D-loop detection

```bash
awk '$4 == "D-loop"' "$ANNO_DIR/annotated_genes.bed"
```

| Observation | Verdict |
|---|---|
| D-loop present with length 500–3,000 bp (typical vertebrate) | ✅ |
| D-loop absent and `--no-trna` was not used | ⚠️ largest intergenic gap < 100 bp — mitogenome is very compact, possibly truncated |
| D-loop length > 5,000 bp | ❌ unusually large — check for missing gene annotations or false positive |
| D-loop on the `-` strand (in the GFF3) | ⚠️ the NCBI GFF3 D-loop uses strand `-` as a convention; not a real issue, but note that `blast_extract_cds.py` always emits strand `+` for the D-loop because it's gap-based. If the reference D-loop has wrap-around coordinates (> genome_size), the GB/GFF parser handles it but you may see start > end in the BED. |

### Check 5 — gene length vs reference

If the reference is available, compare each gene's length in the target to its length in the reference. Length deviations > 20% are biologically suspicious; < 5% is normal for closely related species.

```bash
# Target gene lengths from the FASTA
awk '/^>/{if(name) print name"\t"len; name=$1; sub(/^>/,"",name); len=0; next} {len+=length($1)} END{print name"\t"len}' \
    "$ANNO_DIR/annotated_genes.fasta" > /tmp/target_gene_lengths.tsv

# Compare to user-supplied reference (if available)
# Reference gene lengths can be extracted from --ref_gff/--ref_gb
# Use the same awk pattern against the reference's gene sequences
```

| Per-gene |deviation| | Verdict |
|---|---|
| ≤ 5% (most genes) | ✅ |
| 5–10% (a handful) | ⚠️ biological divergence — typical for congeneric species |
| > 20% (any single gene) | ❌ possible mis-assembly or wrong reference |

### Check 6 — duplicate / overlapping genes

```bash
awk '{print $1"\t"$2"\t"$3}' "$ANNO_DIR/annotated_genes.bed" | sort -k1,1 -k2,2n | \
  awk 'BEGIN{prev_end=-1; prev_name=""} {
    if ($1 == cur_chr && $2 <= prev_end) print "OVERLAP: " prev_name" / "prev_line" and current line";
    cur_chr=$1; prev_end=$3; prev_name=$4; prev_line=$1":"$2"-"$3
  }' "$ANNO_DIR/annotated_genes.bed"
```

(Adjust the script to actually parse the BED file properly.)

| Observation | Verdict |
|---|---|
| No overlaps (allowing ≤ 5 bp for the well-known ATPase8/ATPase6 and ND4L/ND4 overlaps) | ✅ |
| Large overlaps (> 50 bp) between distinct genes | ❌ BLAST hit coordinates overlap — check for paralogous NUMT contamination |

## 3. Output contract — write `annotation-report.md`

Always write to `$RUN_DIR/annotation-report.md`. Template:

```markdown
# BettaMt annotation QC report

Sample:     $SAMPLE_ID
Platform:   $PLATFORM
Pipeline:   $BETA_MT_HOME
Contig:     $POLISHED
Annotation: $ANNO_DIR/annotated_genes.{fasta,gff3,bed}
Generated:  <ISO8601>

## Verdict summary

| Check | Result | Notes |
|---|---|---|
| Gene count | ✅ / ⚠️ / ❌ | <N> genes (<PCGs> PCG, <rRNAs> rRNA, <tRNAs> tRNA, <CR> D-loop) |
| Gene names | ✅ / ⚠️ / ❌ | all canonical |
| Genomic span | ✅ / ⚠️ / ❌ | <annotated_bp> bp |
| D-loop | ✅ / ⚠️ / ❌ | <len> bp at <start>–<end> |
| Length vs reference | ✅ / ⚠️ / ❌ | <per-gene summary> |
| Overlaps | ✅ / ⚠️ / ❌ | <overlaps> |

**Overall**: <PASS / PASS-WITH-WARNINGS / FAIL>

## Gene table

| Gene | Type | Start | End | Strand | Length (bp) | vs ref |
|---|---|---|---|---|---|---|
| 12S | rRNA | … | … | + | … | ±N% |
| … | … | … | … | … | … | … |
| D-loop | CR | … | … | + | … | n/a |

## Notes & recommended follow-ups
- <bullet list of any warnings or failures, with concrete next steps>

## Reproducibility
- params.json: $RUN_DIR/params.json
- pixi.lock:   $BETA_MT_HOME/pixi.lock
- this report: $RUN_DIR/annotation-report.md (regenerable via `bettamt-annotate-qc` skill)
```

## 4. What NOT to do

- Do **not** re-run the pipeline to "fix" annotation issues — produce a report, let the user decide
- Do **not** run `blast_extract_cds.py` directly — that's the pipeline's job
- Do **not** fail the QC on missing tRNAs alone — it could be a taxon model mismatch; suggest re-running with a different `--taxon`
- Do **not** write to `$BETA_MT_HOME/...` paths — only `$RUN_DIR/...`

## 5. Common follow-ups

- If gene count is low → check if the reference annotation is too divergent (use a closer reference or skip annotation)
- If D-loop is missing → check if the polished mitogenome was actually circularized (BettaMt's CIRCLATOR step); if linear, the D-loop is just an intergenic gap and might be < 100 bp
- If tRNAs are missing → re-run with a different `--taxon` (e.g. `mammal` instead of `vertebrate`)
- If gene names are unnormalized → report a bug — the `GENE_ALIASES` dict in `blast_extract_cds.py` should be extended
