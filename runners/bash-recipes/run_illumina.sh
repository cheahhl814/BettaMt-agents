#!/usr/bin/env bash
# run_illumina.sh — BettaMt Illumina mitogenome path, executed manually (bash-recipes runner).
#
# Source of truth: runners/nextflow-runner/subworkflows/illumina_mitogenome.nf
# and modules/local/{qc_short,getorganelle_asm,firstgene,circlator,polca}.nf.
# This script mirrors that wiring command-for-command (see test_smoke.sh):
#
#   QC_SHORT               sequali read-QC report -> fastp trim (adapter detect, front/tail cut, len>=50)
#   GETORGANELLE_ASM       get_organelle_from_reads.py -F animal_mt -R rounds
#   FIRSTGENE              tRNAscan-SE -> awk trnF.bed -> bedtools getfasta
#   CIRCLATOR              circlator fixstart --genes_fa trnF.fasta
#   POLCA                  pypolca run --careful -m 2G  (emits <prefix>_corrected.fasta)
#
# Usage: run_illumina.sh --reads R1 --reads-r2 R2 --ref-mito FASTA [options]
set -euo pipefail

SKILL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_MANIFEST="$SKILL_ROOT/runners/nextflow-runner/pixi.toml"
PIXI() { pixi run --manifest-path "$ENV_MANIFEST" "$@"; }

READS="" READS_R2="" REF="" ROUNDS=15 TAXON="vertebrate" SAMPLE="" OUTDIR="" THREADS=4
while [ $# -gt 0 ]; do case "$1" in
  --reads)     READS="$2";   shift 2;;
  --reads-r2)  READS_R2="$2"; shift 2;;
  --ref-mito)  REF="$2";     shift 2;;
  --rounds)    ROUNDS="$2";  shift 2;;
  --taxon)     TAXON="$2";   shift 2;;
  --sample-id) SAMPLE="$2";  shift 2;;
  --outdir)    OUTDIR="$2";  shift 2;;
  --threads)   THREADS="$2"; shift 2;;
  *) echo "Unknown option: $1"; exit 2;;
esac; done
[ -n "$READS" ] && [ -n "$READS_R2" ] && [ -n "$REF" ] || \
  { echo "ERROR: --reads, --reads-r2 and --ref-mito are required"; exit 2; }
for f in "$READS" "$READS_R2" "$REF"; do [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 2; }; done
SAMPLE="${SAMPLE:-$(basename "$READS" | sed 's/_[Rr]1\.*//; s/\.\(fastq\|fq\)\(\.gz\)\?$//')}"
SAMPLE="${SAMPLE%_*}"
OUTDIR="${OUTDIR:-./bettamt_run_$SAMPLE}"
T="$THREADS"
mkdir -p "$OUTDIR/stages" && cd "$OUTDIR/stages"
echo "[BettaMt Illumina] sample=$SAMPLE threads=$T rounds=$ROUNDS outdir=$(pwd)"

# --- QC_SHORT step 1: sequali read QC -----------------------------------------
echo "[stage 1/5] QC_SHORT (sequali + fastp)"
mkdir -p report
PIXI sequali --threads $T --outdir report \
  "$(readlink -f "$READS")" "$(readlink -f "$READS_R2")"

# --- QC_SHORT step 2: fastp trim ------------------------------------------------
PIXI fastp --in1 "$(readlink -f "$READS")" --in2 "$(readlink -f "$READS_R2")" \
  --out1 trimmed_1.fq.gz --out2 trimmed_2.fq.gz \
  --detect_adapter_for_pe --cut_front --cut_tail --length_required 50 \
  --thread $T --json report/fastp.json --html report/fastp.html

# --- GETORGANELLE_ASM: seed-guided mt assembly ---------------------------------
echo "[stage 2/5] GETORGANELLE_ASM (seed: $(basename "$REF"))"
mkdir -p logs
PIXI get_organelle_from_reads.py \
  -1 "$(readlink -f trimmed_1.fq.gz)" \
  -2 "$(readlink -f trimmed_2.fq.gz)" \
  -s "$(readlink -f "$REF")" \
  -F animal_mt -o "$SAMPLE" -R "$ROUNDS" -t $T \
  --reduce-reads-for-coverage inf --max-reads inf
PIXI sh -c 'cp '"$SAMPLE"'/get_org.log.* logs/ 2>/dev/null || true'
CONTIG="$(ls "$SAMPLE"/*.fasta | head -1)"

# --- FIRSTGENE: trnF extraction -------------------------------------------------
echo "[stage 3/5] FIRSTGENE (taxon model: $TAXON)"
case "$TAXON" in
  mammal)     MODEL_OPT="-M mammal";;
  vertebrate) MODEL_OPT="-M vert";;
  '')         MODEL_OPT="-M vert";;
  *)          MODEL_OPT="-M ${TAXON}";;
esac
PIXI tRNAscan-SE $MODEL_OPT -O -o tRNA.out "$CONTIG"
awk '$5 == "Phe" {
  if ($3 < $4) { print $1"\t"($3-1)"\t"$4"\ttrnF\t"$9"\t+" }
  else         { print $1"\t"($4-1)"\t"$3"\ttrnF\t"$9"\t-" }
}' tRNA.out > trnF.bed
PIXI bedtools getfasta -s -fi "$CONTIG" -bed trnF.bed -fo trnF.fasta

# --- CIRCLATOR: fix start --------------------------------------------------------
echo "[stage 4/5] CIRCLATOR"
PIXI circlator fixstart --genes_fa trnF.fasta "$CONTIG" circlator_fixstart
CIR_CONTIG="circlator_fixstart.fasta"

# --- POLCA: polish with the trimmed short reads -----------------------------------
echo "[stage 5/5] POLCA (pypolca)"
PIXI pypolca run \
  -a "$CIR_CONTIG" \
  -1 trimmed_1.fq.gz \
  -2 trimmed_2.fq.gz \
  -t $T -o pypolca_out -p "$SAMPLE" --careful -m 2G
# pypolca 0.4+ emits <prefix>_corrected.fasta and <prefix>.vcf under -o
# (verified by manual walkthrough 2026-10-18; mirrors runners/nextflow-runner/modules/local/polca.nf).
cp "pypolca_out/${SAMPLE}_corrected.fasta" "../${SAMPLE}_mitogenome.fasta"

echo "[BettaMt Illumina] DONE — polished mitogenome: $OUTDIR/${SAMPLE}_mitogenome.fasta"
echo "Next: qc/bettamt-qc (assembly QC), then optional annotation (annotate_genes.sh)."