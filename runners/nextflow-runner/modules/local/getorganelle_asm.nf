// ===========================================================================
// GETORGANELLE_ASM — assemble mitogenome directly from short (Illumina) reads
//
// GetOrganelle (Kinggerm) was originally written for short-read organelle
// assembly. It performs:
//   1. Seed-based read recruitment with Bowtie2
//   2. Iterative extension of organelle reads
//   3. SPAdes de-Bruijn assembly at multiple k-mer sizes
//   4. Graph disentanglement and target contig selection
//
// This is the *primary assembler* in the short-read path. The output FASTA
// (typically <output>.fasta) is the assembled mitogenome.
//
// Required parameters:
//   - seed    : reference mitogenome (FASTA) — used as initial seed
//   - R       : max extension rounds (10-15 typical for animal_mt)
// ===========================================================================
process GETORGANELLE_ASM {
    tag        "GetOrganelle: Assemble mitogenome from short reads"
    publishDir "${launchDir}/results", mode: 'copy', overwrite: false, pattern: '**'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path  seed
    val   rounds

    output:
    tuple val(sample_id), path("${sample_id}/*.fasta"), emit: contig
    path "logs",                                         emit: logs

    script:
    def PIXI = "pixi run --manifest-path ${baseDir}/pixi.toml"
    // --reduce-reads-for-coverage inf + --max-reads inf: avoid GetOrganelle's
    // default 200x cap, which can be too aggressive for deeply sequenced
    // whole-genome libraries. See GetOrganelle issue #116.
    """
    set -euo pipefail
    mkdir -p logs

    ${PIXI} get_organelle_from_reads.py \\
        -1 ${r1} \\
        -2 ${r2} \\
        -s ${seed} \\
        -F animal_mt \\
        -o ${sample_id} \\
        -R ${rounds} \\
        -t ${task.cpus} \\
        --reduce-reads-for-coverage inf \\
        --max-reads               inf

    # Surface the run log for debugging even when assembly succeeds
    cp ${sample_id}/get_org.log.* logs/ 2>/dev/null || true
    """
}
