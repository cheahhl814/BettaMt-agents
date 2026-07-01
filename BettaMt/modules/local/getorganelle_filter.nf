// ===========================================================================
// GETORGANELLE_FILTER — extract the mitochondrial contig from an assembly graph
//
// Used by the long-read (ONT) path after Flye. GetOrganelle is run as an
// assembly-graph filter here, not a full assembler.
// ===========================================================================
process GETORGANELLE_FILTER {
    tag        "Filter mitochondrial contigs from assembly graph"
    publishDir "${launchDir}/results", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path gfa

    output:
    path "circlator_clean.fasta", emit: contig

    script:
    """
    pixi run --manifest-path ${baseDir}/pixi.toml get_organelle_config --add animal_mt
    pixi run --manifest-path ${baseDir}/pixi.toml get_organelle_from_assembly \\
        -F animal_mt -g ${gfa} -o get_organelle -t ${task.cpus}
    get_org_fasta=\$(ls get_organelle/*.fasta)
    pixi run --manifest-path ${baseDir}/pixi.toml \\
        circlator clean --min_contig_length 12000 \${get_org_fasta} circlator_clean
    """
}
