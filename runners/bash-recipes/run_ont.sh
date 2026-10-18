#!/usr/bin/env bash
# run_ont.sh — BettaMt ONT mitogenome path, executed manually (bash-recipes runner).
#
# Source of truth: runners/nextflow-runner/subworkflows/ont_mitogenome.nf and
# runners/nextflow-runner/modules/local/*.nf. This script mirrors that wiring
# command-for-command (see test_smoke.sh for the enforced lock-step check):
#
#   BAIT_MITO              seqkit seq filter -> bbduk bait -> minimap2 map-ont -> samtools remap->fastq
#   FLYE                   flye --meta --nano-raw
#   GETORGANELLE_FILTER    get_organelle_config --add animal_mt -> get_organelle_from_assembly -F animal_mt -> circlator clean 12000
#   FIRSTGENE              tRNAscan-SE (-M model) -> awk trnF.bed -> bedtools getfasta
#   CIRCLATOR              circlator fixstart --genes_fa trnF.fasta
#   RACON                  minimap2 map-ont -> racon -> header rename
#
# All tool invocations go through pixi against runners/nextflow-runner/pixi.toml
# so versions stay pinned to the same environment the Nextflow runner uses.
#
# Usage: run_ont.sh --reads FASTQ --ref-mito FASTA [options]
set -euo pipefail

SKILL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_MANIFEST="$SKILL_ROOT/runners/nextflow-runner/pixi.toml"
PIXI() { pixi run --manifest-path "$ENV_MANIFEST" "$@"; }

READS="" REF="" SIZE="16k" TAXON="vertebrate" SAMPLE="" OUTDIR="" THREADS=4
MIN_LEN=14000 MAX_LEN=20000 SKIP_TRNA=0
while [ $# -gt 0 ]; do case "$1" in
  --reads)     READS="$2";   shift 2;;
  --ref-mito)  REF="$2";     shift 2;;
  --size)      SIZE="$2";    shift 2;;
  --taxon)     TAXON="$2";   shift 2;;
  --sample-id) SAMPLE="$2";  shift 2;;
  --outdir)    OUTDIR="$2";  shift 2;;
  --threads)   THREADS="$2"; shift 2;;
  --skip-annotation-prep) SKIP_TRNA=1; shift;;   # reserved: skip tRNA model prep only
  *) echo "Unknown option: $1"; exit 2;;
esac; done
[ -n "$READS" ] && [ -n "$REF" ] || { echo "ERROR: --reads and --ref-mito are required"; exit 2; }
[ -f "$READS" ] || { echo "ERROR: reads not found: $READS";  exit 2; }
[ -f "$REF"   ] || { echo "ERROR: ref not found: $REF";      exit 2; }
SAMPLE="${SAMPLE:-$(basename "$READS" | sed 's/\.\(fastq\|fq\)\(\.gz\)\?$//')}"
OUTDIR="${OUTDIR:-./bettamt_run_$SAMPLE}"
T="$THREADS"
mkdir -p "$OUTDIR/stages" && cd "$OUTDIR/stages"
echo "[BettaMt ONT] sample=$SAMPLE threads=$T outdir=$(pwd)"

# --- BAIT_MITO: length filter -> bait -> remap -> fastq ----------------------
echo "[stage 1/6] BAIT_MITO"
PIXI seqkit seq -m $MIN_LEN -M $MAX_LEN "$READS" > "${SAMPLE}_filtered.fastq"
PIXI bbduk.sh in="${SAMPLE}_filtered.fastq" outm="${SAMPLE}_baited.fastq" \
    ref="$(readlink -f "$REF")" k=31 hdist=1 threads=$T
PIXI minimap2 -ax map-ont -t $T "$(readlink -f "$REF")" "${SAMPLE}_baited.fastq" | \
  PIXI samtools view -bF 4 - | PIXI samtools fastq - > "${SAMPLE}_mapped.fastq"

# --- FLYE: assemble ----------------------------------------------------------
echo "[stage 2/6] FLYE"
PIXI flye --genome-size "$SIZE" --meta --threads $T --out-dir flye --nano-raw "${SAMPLE}_mapped.fastq"

# --- GETORGANELLE_FILTER: mt contigs from the graph -> circlator clean -------
echo "[stage 3/6] GETORGANELLE_FILTER"
# NOTE: entry points are .py-suffixed in the pinned env (get_organelle_config.py,
# get_organelle_from_assembly.py) — runners/nextflow-runner/modules/local/getorganelle_filter.nf
# calls them WITHOUT .py and will fail against this env (latent NF bug; use this recipe).
PIXI get_organelle_config.py --add animal_mt
PIXI get_organelle_from_assembly.py -F animal_mt -g flye/assembly_graph.gfa -o get_organelle -t $T
GET_ORG_FASTA="$(ls get_organelle/*.fasta)"
PIXI circlator clean --min_contig_length 12000 "$GET_ORG_FASTA" circlator_clean

# --- FIRSTGENE: trnF extraction (drives fixstart) ----------------------------
echo "[stage 4/6] FIRSTGENE (taxon model: $TAXON)"
case "$TAXON" in
  mammal)     MODEL_OPT="-M mammal";;
  vertebrate) MODEL_OPT="-M vert";;
  '')         MODEL_OPT="-M vert";;
  *)          MODEL_OPT="-M ${TAXON}";;
esac
PIXI tRNAscan-SE $MODEL_OPT -O -o tRNA.out circlator_clean.fasta
awk '$5 == "Phe" {
  if ($3 < $4) { print $1"\t"($3-1)"\t"$4"\ttrnF\t"$9"\t+" }
  else         { print $1"\t"($4-1)"\t"$3"\ttrnF\t"$9"\t-" }
}' tRNA.out > trnF.bed
PIXI bedtools getfasta -s -fi circlator_clean.fasta -bed trnF.bed -fo trnF.fasta

# --- CIRCLATOR: fix start against trnF ---------------------------------------
echo "[stage 5/6] CIRCLATOR"
PIXI circlator fixstart --genes_fa trnF.fasta circlator_clean.fasta circlator_fixstart

# --- RACON: polish with the baited long reads --------------------------------
echo "[stage 6/6] RACON"
mv circlator_fixstart.fasta circular.fasta
PIXI minimap2 -ax map-ont -t $T circular.fasta "${SAMPLE}_mapped.fastq" > "${SAMPLE}.sam"
PIXI racon -t $T "${SAMPLE}_mapped.fastq" "${SAMPLE}.sam" circular.fasta > "${SAMPLE}_racon_tmp.fasta"
sed "s/^>.*/>${SAMPLE}_mitochondrion/" "${SAMPLE}_racon_tmp.fasta" > "${SAMPLE}_racon.fasta"

ln -sfn "$(readlink -f "${SAMPLE}_racon.fasta")" "$OUTDIR/${SAMPLE}_mitogenome.fasta"
echo "[BettaMt ONT] DONE — polished mitogenome: $OUTDIR/${SAMPLE}_mitogenome.fasta"
echo "Next: qc/bettamt-qc (assembly QC), then optional annotation (annotate_genes.sh)."