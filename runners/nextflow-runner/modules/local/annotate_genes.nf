// ===========================================================================
// ANNOTATE_GENES — homology-based annotation of all mitochondrial genes
//
// Wraps blast_extract_cds.py to produce:
//   * multi-gene FASTA     (per-gene sequences, sorted by position)
//   * GFF3 annotation      (standard gene annotations)
//   * BED coordinates      (genomic intervals)
//
// Uses a reference annotation (GFF3 + FASTA, or GenBank) to define query
// sequences for BLASTN, then transfers annotations to the target mitogenome
// coordinates.  Also runs tRNAscan-SE for tRNA genes and detects the D-loop
// (control region) as the largest intergenic gap.
//
// Input:
//   polished      — final circularized, rotated mitogenome FASTA
//   ref_gff + ref_fasta  — GFF3 reference annotation + reference FASTA, OR
//   ref_gb        — GenBank reference annotation (self-contained)
//   trnascan_path — optional path to tRNAscan-SE binary (auto-detected otherwise)
//   skip_trna     — if true, skip tRNA annotation (PCGs + rRNAs only)
//
// Output:
//   .fasta  — multi-gene FASTA
//   .gff3   — gene annotations
//   .bed    — genomic coordinates
// ===========================================================================
process ANNOTATE_GENES {
    tag "Annotate mitochondrial genes (BLAST + tRNAscan-SE)"
    publishDir "${launchDir}/results/annotation", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path polished
    path ref_gff
    path ref_fasta
    path ref_gb
    path trnascan_path
    val  skip_trna

    output:
    path "annotated_genes.fasta", emit: fasta
    path "annotated_genes.gff3",  emit: gff3
    path "annotated_genes.bed",   emit: bed

    script:
    def PIXI     = "pixi run --manifest-path ${baseDir}/pixi.toml"
    def SCRIPT   = "${baseDir}/../bin/blast_extract_cds.py"
    def args     = "--target ${polished} --out annotated_genes"
    def ref_args = ""
    if (ref_gb) {
        ref_args = "--ref-gb ${ref_gb}"
    } else if (ref_gff && ref_fasta) {
        ref_args = "--ref-gff ${ref_gff} --ref-fasta ${ref_fasta}"
    } else if (ref_gff) {
        ref_args = "--ref-gff ${ref_gff}"
    }
    def trna_args = skip_trna ? "--no-trna" : (trnascan_path ? "--trnascan ${trnascan_path}" : "")
    """
    ${PIXI} python ${SCRIPT} \\
        ${args} ${ref_args} ${trna_args}
    """
}
