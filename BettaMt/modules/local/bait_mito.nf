// ===========================================================================
// BAIT_MITO — filter ONT reads by length and bait against a reference mitogenome
//
// Inputs : FASTQ (ONT), reference mitochondrial FASTA
// Output : FASTQ of reads that mapped to the reference (mt-enriched)
//
// Strategy:
//   1. seqkit filters reads to a 14–20 kb window (vertebrate mitogenome size).
//   2. bbduk k-mer matches the windowed reads against the reference (outm=k-mer hits).
//   3. minimap2 (map-ont) maps those k-mer hits back; -F 4 keeps only mapped reads.
//   4. samtools fastq emits the mt-enriched FASTQ.
// ===========================================================================
process BAIT_MITO {
    tag        "Baiting for mitochondrial reads"
    publishDir "${launchDir}/results/bait", mode: 'copy', overwrite: false, pattern: '**'

    input:
    path fastq
    path mt_ref

    output:
    path "${fastq.baseName}_mapped.fastq", emit: baited

    script:
    def PIXI = "pixi run --manifest-path ${baseDir}/pixi.toml"
    """
    ${PIXI} seqkit seq -m 14000 -M 20000 ${fastq} > ${fastq.baseName}_filtered.fastq
    ${PIXI} bbduk in=${fastq.baseName}_filtered.fastq outm=${fastq.baseName}_baited.fastq ref=${mt_ref} k=31 hdist=1 threads=${task.cpus}
    ${PIXI} minimap2 -ax map-ont -t ${task.cpus} ${mt_ref} ${fastq.baseName}_baited.fastq | \\
        ${PIXI} samtools view -bF 4 - | \\
        ${PIXI} samtools fastq - > ${fastq.baseName}_mapped.fastq
    """
}
