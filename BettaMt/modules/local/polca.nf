// ===========================================================================
// POLCA — short-read polishing of the final circularized contig
//
// pypolca (gbouras13) is a Python re-implementation of POLCA from MaSuRCA
// (Zimin & Salzberg, PLoS Comput Biol 2020). It:
//   1. Aligns Illumina reads to the assembly with BWA-MEM
//   2. Sorts/indexes with samtools
//   3. Calls variants with freebayes
//   4. Applies consensus corrections above configurable thresholds
//
// We pass --careful (min_alt 4, min_ratio 3) — recommended for bacterial
// isolate polishing and an appropriate balance for organelle genomes.
// ===========================================================================
process POLCA {
    tag        "POLCA: Polish mitogenome with ${sample_id} Illumina reads"
    publishDir "${launchDir}/results/polish", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(sample_id), path(contig), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_polca.fasta"), emit: polished
    path "pypolca_out",                                     emit: vcf

    script:
    def PIXI = "pixi run --manifest-path ${baseDir}/pixi.toml"
    """
    set -euo pipefail

    # pypolca emits ${prefix}.fasta and ${prefix}.vcf under -o.
    ${PIXI} pypolca run \\
        -a ${contig} \\
        -1 ${r1} \\
        -2 ${r2} \\
        -t ${task.cpus} \\
        -o pypolca_out \\
        -p ${sample_id} \\
        --careful \\
        -m 2G

    cp pypolca_out/${sample_id}.fasta ${sample_id}_polca.fasta
    """
}
