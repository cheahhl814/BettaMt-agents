// ===========================================================================
// FLYE — long-read de novo assembly (ONT mode)
//
// Inputs : mt-enriched FASTQ, estimated genome size (e.g. "16k")
// Outputs: assembly.fasta and assembly_graph.gfa
// ===========================================================================
process FLYE {
    tag        "Flye: Assemble the mitochondrial genome"
    publishDir "${launchDir}/results", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path fastq
    val size

    output:
    path "flye"
    path "flye/assembly.fasta",                emit: contig
    path "flye/assembly_graph.gfa",            emit: gfa

    script:
    """
    pixi run --manifest-path ${baseDir}/pixi.toml \\
        flye --genome-size ${size} --meta --threads ${task.cpus} --out-dir flye --nano-raw ${fastq}
    """
}
