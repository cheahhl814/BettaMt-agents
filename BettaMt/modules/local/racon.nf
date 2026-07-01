// ===========================================================================
// RACON — long-read polishing of the circularized contig
//
// Used only by the ONT subworkflow. Races minimap2 (map-ont) alignments of
// raw ONT reads back to the rotated contig and emits a polished FASTA.
// ===========================================================================
process RACON {
    tag        "Racon: Polish mitochondrial genome (long reads)"
    publishDir "${launchDir}/results/polish", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path contig
    path fastq

    output:
    path "*_racon.fasta"

    script:
    def sample_id = fastq.baseName
    def PIXI      = "pixi run --manifest-path ${baseDir}/pixi.toml"
    """
    mv ${contig} circular.fasta
    ${PIXI} minimap2 -ax map-ont -t ${task.cpus} circular.fasta ${fastq} > ${sample_id}.sam
    ${PIXI} racon -t ${task.cpus} ${fastq} ${sample_id}.sam circular.fasta > ${sample_id}_racon_tmp.fasta
    sed "s/^>.*/>${sample_id}_mitochondrion/" ${sample_id}_racon_tmp.fasta > ${sample_id}_racon.fasta
    """
}
