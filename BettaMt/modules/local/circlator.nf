// ===========================================================================
// CIRCLATOR — circularize and rotate the contig to a canonical start site
//
// Shared by both ONT and Illumina subworkflows. Reads "circlator_clean.fasta"
// written by the upstream filter/assembly step in the same work directory.
// ===========================================================================
process CIRCLATOR {
    tag        "Circlator: Circularization and Rotation"
    publishDir "${launchDir}/results/circlator", mode: 'copy', overwrite: false

    input:
    path contig
    path genes_fa

    output:
    path "circlator_fixstart.fasta", emit: contig

    script:
    """
    pixi run --manifest-path ${baseDir}/pixi.toml \\
        circlator fixstart --genes_fa ${genes_fa} circlator_clean.fasta circlator_fixstart
    """
}
