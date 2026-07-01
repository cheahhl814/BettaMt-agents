---
name: bettamt-qc
description: Quality-check a finished BettaMt mitogenome assembly. Use after a successful nextflow run when the user wants to know if the assembly is complete, circular, NUMT-free, and publishable. Triggers: "bettamt QC", "check my mitogenome", "is the mitogenome complete", "post-assembly QC", "bettamt report", "did the assembly work".
---

# bettamt-qc — BettaMt post-assembly quality control

This skill runs **after** a successful `nextflow run`. It never re-runs the pipeline. It produces one `report.md` with verdicts, evidence, and (optionally) a coverage profile plot.

## 0. Locate the pipeline and the run

```bash
: "${BETA_MT_HOME:=$(cd "$(dirname "$SKILL_DIR")/../BettaMt" 2>/dev/null && pwd)}"
test -f "$BETA_MT_HOME/main.nf" || { echo "Cannot locate BettaMt. Set BETA_MT_HOME."; exit 1; }

# RUN_DIR is wherever the user launched nextflow
RUN_DIR="${RUN_DIR:-$BETA_MT_HOME}"
test -d "$RUN_DIR/results" || { echo "No results/ directory in $RUN_DIR"; exit 1; }
```

## 1. Identify platform and final contig

```bash
PARAMS="$RUN_DIR/params.json"
PLATFORM=$(jq -r '.platform' < "$PARAMS")

case "$PLATFORM" in
  ont)
    POLISHED=$(ls "$RUN_DIR"/results/polish/*_racon.fasta 2>/dev/null | head -1)
    TRIMMED_OR_BAITED="$RUN_DIR/results/bait"
    ;;
  illumina)
    POLISHED=$(ls "$RUN_DIR"/results/polish/*_polca.fasta 2>/dev/null | head -1)
    TRIMMED_OR_BAITED="$RUN_DIR/results/qc"  # trim reports
    ;;
  *) echo "Unknown platform: $PLATFORM"; exit 1 ;;
esac

test -n "$POLISHED" || { echo "Polished FASTA not found under results/polish/"; exit 1; }
```

## 2. Run the QC checks

Each check produces a verdict (✅ pass, ⚠️ warn, ❌ fail) and a one-line reason.

### Check 1 — contig size and contiguity

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" seqkit fx2tab -n -l "$POLISHED"
```

| Length | Verdict |
|---|---|
| 14 000–20 000 bp (typical vertebrate mt) | ✅ |
| 10 000–14 000 or 20 000–25 000 | ⚠️ possible duplication or truncation |
| < 10 000 or > 25 000 | ❌ likely not a mitogenome (or NUMT) |
| Multiple contigs of similar length | ⚠️ fragmented — verify with `grep -c '>' "$POLISHED"` |

### Check 2 — circularity

```bash
# GetOrganelle outputs typically have a 75-bp terminal repeat (Hamming distance ~0)
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" seqkit head -n 2 "$POLISHED"  # first 2 records
# Compare first 75 bp of record 1 vs last 75 bp of last record
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" seqkit subseq -r 1:75 "$POLISHED"
# A circular mitogenome should show the same 75 bp at start and end
```

| Observation | Verdict |
|---|---|
| First 75 bp == last 75 bp (exact or near-exact) | ✅ circularized |
| First 75 bp ≈ last 75 bp with 1–3 mismatches | ✅ circular (Racon/POLCA polish artifact) |
| No obvious repeat | ⚠️ may be linear — check upstream contig |

### Check 3 — tRNA gene set

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" tRNAscan-SE -M vert -O -o /tmp/qc-trna.out "$POLISHED"
awk '$5 != "" {print $5}' /tmp/qc-trna.out | sort -u
```

| tRNAs found | Verdict |
|---|---|
| 22 unique tRNAs (vertebrate standard) | ✅ complete |
| 20–21 | ⚠️ nearly complete; missing tRNA could be a real biology or a model miss |
| < 18 | ❌ incomplete — likely a fragment or wrong model |
| Run with `-M anonym` if < 15 found | re-check; this should always pass after re-run |

### Check 4 — coverage profile (NUMT risk)

Map a downsampled subset of reads back to the polished contig and look for sudden coverage jumps (NUMT signature = a 5–10× local spike):

```bash
SAMPLE_READS=""
if [ "$PLATFORM" = "illumina" ]; then
    SAMPLE_READS=$(ls "$RUN_DIR"/results/qc/trimmed_1.fq.gz 2>/dev/null)
else
    SAMPLE_READS=$(ls "$RUN_DIR"/results/bait/*_mapped.fastq 2>/dev/null)
fi

pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" minimap2 -ax $([ "$PLATFORM" = "illumina" ] && echo sr || echo map-ont) -t 4 "$POLISHED" "$SAMPLE_READS" \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools view -bF 4 - \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools sort -o /tmp/qc.bam -
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools index /tmp/qc.bam
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools depth -a /tmp/qc.bam > /tmp/qc.depth

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

| Observation | Verdict |
|---|---|
| max coverage ≤ 5× median | ✅ clean |
| 5–10× median, single window | ⚠️ possible NUMT — inspect with IGV |
| > 10× median, single window | ❌ strong NUMT signal — investigate that locus |

If the user wants a plot, generate a simple ASCII histogram of the depth file (or a PNG with `python -c "import matplotlib..."` if `matplotlib` is available — it is **not** in the BettaMt pixi env by default, so prefer ASCII).

### Check 5 — identity to reference

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" minimap2 -ax $([ "$PLATFORM" = "illumina" ] && echo sr || echo map-ont) "$REF_MITO" "$POLISHED" \
  | pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" samtools view - \
  | awk '$6 != "*" { n++; m+=length($10)*($1/$1); } END { print "ref_pct_identity=" m/(n*length($10)) }'
# Simpler: just report alignment length and number of substitutions from CIGAR
```

| Identity to ref | Verdict |
|---|---|
| > 95% | ✅ expected (same / closely related species) |
| 85–95% | ⚠️ congeneric — verify expected divergence |
| < 85% | ❌ suspicious — might be a different locus |

### Check 6 — GC content and skew

```bash
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" seqkit fx2tab -n -l -g "$POLISHED"
```

| GC% | Verdict |
|---|---|
| 35–50% (most vertebrates) | ✅ |
| 25–35 or 50–60% | ⚠️ unusual but possible (e.g. AT-rich fish mitogenomes) |
| < 25% or > 60% | ❌ suspicious |

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
pixi run --manifest-path "$BETA_MT_HOME/pixi.toml" python -c "
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
