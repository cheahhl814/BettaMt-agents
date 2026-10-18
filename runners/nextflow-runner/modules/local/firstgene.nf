// ===========================================================================
// FIRSTGENE — annotate trnF (Phe) and extract its sequence for rotation
//
// Shared by both ONT and Illumina subworkflows. The output trnF.fasta is
// consumed by CIRCLATOR fixstart to choose a canonical start position.
// ===========================================================================
process FIRSTGENE {
    tag "Initial annotation and extraction of the first gene (trnF)"

    input:
    path contig
    val taxon

    output:
    path "trnF.fasta", emit: fasta

    script:
    // Pick tRNAscan-SE covariance model. Default to vertebrate when unset.
    def model_opt = ""
    if (taxon == 'mammal') {
        model_opt = "-M mammal"
    } else if (taxon == 'vertebrate' || taxon == '') {
        model_opt = "-M vert"
    } else {
        // Pass through any other model name the user supplies
        model_opt = "-M ${taxon}"
    }

    def PIXI = "pixi run --manifest-path ${baseDir}/pixi.toml"
    """
    ${PIXI} tRNAscan-SE ${model_opt} -O -o tRNA.out ${contig}
    awk '\$5 == "Phe" {
        if (\$3 < \$4) { print \$1"\\t"(\$3-1)"\\t"\$4"\\ttrnF\\t"\$9"\\t+" }
        else          { print \$1"\\t"(\$4-1)"\\t"\$3"\\ttrnF\\t"\$9"\\t-" }
    }' tRNA.out > trnF.bed
    ${PIXI} bedtools getfasta -s -fi ${contig} -bed trnF.bed -fo trnF.fasta
    """
}
