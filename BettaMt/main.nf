#!/usr/bin/env nextflow

// ===========================================================================
// BettaMt — entry point
//
// Two interchangeable subworkflows, selected by --platform:
//   --platform ont       Oxford Nanopore long reads
//   --platform illumina  Paired-end short reads (Illumina / Element / MGI)
//
// Shared parameters:
//   --taxon     tRNAscan-SE model: mammal | vertebrate | (any -M name)
//   --size      Estimated mitogenome size, ONT path only (e.g. "16k")
//   --ref_mito  Reference mitogenome (ONT: bait;  Illumina: assembly seed)
//
// Platform-specific parameters:
//   --reads       ONT: single FASTQ        | Illumina: R1 FASTQ
//   --reads_r2    Illumina only: R2 FASTQ
//   --sample_id   Illumina only: sample identifier (defaults to reads basename)
//   --rounds      Illumina only: GetOrganelle -R value (default 15)
// ===========================================================================

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Parameters
// -----------------------------------------------------------------------------
params.platform  = 'ont'                       // 'ont' | 'illumina'
params.reads     = ''                          // FASTQ (ONT) or R1 (Illumina)
params.reads_r2  = ''                          // R2 (Illumina only)
params.ref_mito  = ''                          // seed/bait reference mitogenome
params.size      = '16k'                       // ONT only: expected mitogenome size
params.taxon     = 'vertebrate'                // tRNAscan-SE model
params.sample_id = null                        // Illumina only: sample name
params.rounds    = 15                          // Illumina only: GetOrganelle -R

// -----------------------------------------------------------------------------
// Subworkflows
// -----------------------------------------------------------------------------
include { ONT_MITOGENOME      } from './subworkflows/ont_mitogenome.nf'
include { ILLUMINA_MITOGENOME } from './subworkflows/illumina_mitogenome.nf'

// -----------------------------------------------------------------------------
// Main workflow
// -----------------------------------------------------------------------------
workflow {

    if (params.platform == 'ont') {
        if (!params.reads)   { error "ERROR: --reads (ONT FASTQ) is required for --platform ont" }
        if (!params.ref_mito){ error "ERROR: --ref_mito is required for --platform ont" }

        ch_reads    = channel.fromPath(params.reads,    checkIfExists: true)
        ch_ref_mito = channel.fromPath(params.ref_mito, checkIfExists: true)
        size        = channel.value(params.size)
        taxon       = channel.value(params.taxon)

        ONT_MITOGENOME(ch_reads, ch_ref_mito, size, taxon)

    } else if (params.platform == 'illumina') {
        if (!params.reads)    { error "ERROR: --reads (R1 FASTQ) is required for --platform illumina" }
        if (!params.reads_r2) { error "ERROR: --reads_r2 (R2 FASTQ) is required for --platform illumina" }
        if (!params.ref_mito) { error "ERROR: --ref_mito is required as a seed for --platform illumina" }

        def sid    = params.sample_id ?: file(params.reads).baseName
                            .replaceAll('_(R1|1)\\.(fastq|fq)(\\.gz)?$', '')
        ch_r1      = channel.fromPath(params.reads,    checkIfExists: true)
        ch_r2      = channel.fromPath(params.reads_r2, checkIfExists: true)
        ch_paired  = channel.value(sid)
                            .combine(ch_r1)
                            .combine(ch_r2)
        ch_seed    = channel.fromPath(params.ref_mito, checkIfExists: true)
        rounds     = channel.value(params.rounds)
        taxon      = channel.value(params.taxon)

        ILLUMINA_MITOGENOME(ch_paired, ch_seed, rounds, taxon)

    } else {
        error "ERROR: --platform must be 'ont' or 'illumina' (got '${params.platform}')"
    }
}

workflow.onComplete {
    log.info "BettaMt finished. Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Results: ${launchDir}/results"
}
