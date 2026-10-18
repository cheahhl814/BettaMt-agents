#!/usr/bin/env bash
# annotate_genes.sh — optional post-polish gene annotation (BettaMt bash-recipes runner).
#
# Mirrors runners/nextflow-runner/modules/local/annotate_genes.nf command-for-command:
# it wraps runners/nextflow-runner/bin/blast_extract_cds.py through the pinned pixi env.
#
# Run AFTER a successful run_ont.sh / run_illumina.sh, pointing --target at the
# polished mitogenome. Requires a reference annotation: GFF3 (optionally +
# companion FASTA) or a GenBank file — exactly as --ref_gff/--ref_gff_fasta/--ref_gb
# do for the pipeline.
#
# Usage: annotate_genes.sh --target Mitogenome.fasta --ref-gff ref.gff3 [--ref-fasta ref.fasta]
#        annotate_genes.sh --target Mitogenome.fasta --ref-gb ref.gb [--no-trna | --trnascan PATH]
set -euo pipefail

SKILL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_MANIFEST="$SKILL_ROOT/runners/nextflow-runner/pixi.toml"
PIXI() { pixi run --manifest-path "$ENV_MANIFEST" "$@"; }
SCRIPT="$SKILL_ROOT/runners/nextflow-runner/bin/blast_extract_cds.py"

TARGET="" REF_GFF="" REF_FASTA="" REF_GB="" TRNASCAN="" NO_TRNA=0 OUT="annotated_genes"
while [ $# -gt 0 ]; do case "$1" in
  --target)    TARGET="$2";    shift 2;;
  --ref-gff)   REF_GFF="$2";   shift 2;;
  --ref-fasta) REF_FASTA="$2"; shift 2;;
  --ref-gb)    REF_GB="$2";    shift 2;;
  --trnascan)  TRNASCAN="$2";  shift 2;;
  --no-trna)   NO_TRNA=1;      shift;;
  --out)       OUT="$2";       shift 2;;
  *) echo "Unknown option: $1"; exit 2;;
esac; done
[ -n "$TARGET" ] && [ -f "$TARGET" ] || { echo "ERROR: --target (polished mitogenome fasta) required"; exit 2; }
[ -n "$REF_GFF" ] || [ -n "$REF_GB" ] || { echo "ERROR: provide --ref-gff (optionally + --ref-fasta) or --ref-gb"; exit 2; }

ARGS=(--target "$TARGET" --out "$OUT")
if [ -n "$REF_GB" ]; then
  ARGS+=(--ref-gb "$REF_GB")
else
  ARGS+=(--ref-gff "$REF_GFF")
  [ -n "$REF_FASTA" ] && ARGS+=(--ref-fasta "$REF_FASTA")
fi
if [ "$NO_TRNA" -eq 1 ]; then ARGS+=(--no-trna)
elif [ -n "$TRNASCAN" ]; then ARGS+=(--trnascan "$TRNASCAN"); fi

echo "[annotate] $0 on $(basename "$TARGET") -> $OUT/"
PIXI python "$SCRIPT" "${ARGS[@]}"
echo "[annotate] outputs: ${OUT}.fasta ${OUT}.gff3 ${OUT}.bed (verify with qc/bettamt-qc annotation checks)"