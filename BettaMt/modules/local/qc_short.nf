// ===========================================================================
// QC_SHORT — quality control + adapter/quality trimming for short (Illumina) reads
//
// Strategy:
//   - sequali  -> per-base & per-read QC report (HTML + text summary)
//   - fastp    -> adapter trimming, quality filtering, length filtering
//                 (default keeps reads >= 50 bp after trimming)
//
// Output:
//   - report/        : sequali HTML + summary  (publishDir)
//   - trimmed_1.fq.gz, trimmed_2.fq.gz : fastp output fed into the assembler
// ===========================================================================
process QC_SHORT {
    tag        "QC + trim: ${sample_id}"
    publishDir "${launchDir}/results/qc", mode: 'copy', overwrite: false, pattern: 'report/*'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("trimmed_1.fq.gz"), path("trimmed_2.fq.gz"), emit: trimmed
    path "report",                                                          emit: report

    script:
    def PIXI = "pixi run --manifest-path ${baseDir}/pixi.toml"
    """
    set -euo pipefail
    mkdir -p report

    # 1. Read QC (sequali) — produces per-base plots and a tabular summary
    ${PIXI} sequali \\
        --threads ${task.cpus} \\
        --outdir  report \\
        --force \\
        ${r1} ${r2}

    # 2. Adapter + quality + length trimming (fastp)
    #    - detect_adapters lets fastp infer Illumina/TruSeq/ Nextera adapters
    #    - cut_front / cut_tail trim each end regardless of quality
    #    - length_required 50 drops fragments too short to span a k-mer
    ${PIXI} fastp \\
        --in1           ${r1} \\
        --in2           ${r2} \\
        --out1          trimmed_1.fq.gz \\
        --out2          trimmed_2.fq.gz \\
        --detect_adapter_for_pe \\
        --cut_front --cut_tail \\
        --length_required 50 \\
        --thread        ${task.cpus} \\
        --json          report/fastp.json \\
        --html          report/fastp.html
    """
}
