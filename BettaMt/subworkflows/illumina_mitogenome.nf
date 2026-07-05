// ===========================================================================
// ILLUMINA_MITOGENOME — subworkflow for short-read (Illumina) assembly
//
//   (r1,r2)       ─► QC_SHORT              ─► trimmed
//   trimmed, seed, rounds ─► GETORGANELLE_ASM ─► assembled
//   assembled           ─► FIRSTGENE        ─► trnF
//   assembled, trnF     ─► CIRCLATOR        ─► rotated
//   rotated, trimmed    ─► POLCA            ─► polished (final mitogenome)
//
// Notes:
//   * --seed is the *initial seed* FASTA (typically a related-species
//     mitogenome). It is required because pure de novo organelle assembly
//     from a whole-genome library is rarely successful.
//   * --rounds sets GetOrganelle's maximum extension iterations. 10–15 is
//     a safe default for animal_mt.
// ===========================================================================
nextflow.enable.dsl = 2

include { QC_SHORT          } from '../modules/local/qc_short.nf'
include { GETORGANELLE_ASM } from '../modules/local/getorganelle_asm.nf'
include { FIRSTGENE        } from '../modules/local/firstgene.nf'
include { CIRCLATOR        } from '../modules/local/circlator.nf'
include { POLCA            } from '../modules/local/polca.nf'
include { ANNOTATE_GENES   } from '../modules/local/annotate_genes.nf'

workflow ILLUMINA_MITOGENOME {

    take:
    ch_reads_paired  // channel: [ val(sample_id), path(r1), path(r2) ]
    ch_seed          // channel: [ path(seed_fasta) ]
    ch_ref_gff       // channel: [ path(gff3) ]  — optional, for gene annotation
    ch_ref_fasta     // channel: [ path(fasta) ] — optional, companion to ref_gff
    ch_ref_gb        // channel: [ path(gb) ]    — optional, alternative to gff+fasta
    ch_trnascan      // channel: [ path(binary) ] — optional
    skip_trna        // value:    boolean
    rounds           // value:    GetOrganelle -R value
    taxon            // value:    tRNAscan-SE model (e.g. "vertebrate")

    main:
    QC_SHORT(ch_reads_paired)
    GETORGANELLE_ASM(QC_SHORT.out.trimmed, ch_seed, rounds)
    FIRSTGENE(GETORGANELLE_ASM.out.contig, taxon)

    // CIRCLATOR expects a single contig path; flatten the [id, fasta] tuple
    ch_clean = GETORGANELLE_ASM.out.contig.map { sid, fa -> fa }
    CIRCLATOR(ch_clean, FIRSTGENE.out.fasta)

    // POLCA joins the rotated contig with the trimmed reads on sample_id
    ch_polca_in = CIRCLATOR.out.contig
        .map { fa -> [ fa.baseName.replaceAll('_circlator_fixstart$',''), fa ] }
        .join(QC_SHORT.out.trimmed, by: 0)
    POLCA(ch_polca_in)

    // Gene annotation — only if a reference annotation was provided
    if (ch_ref_gff || ch_ref_gb) {
        ANNOTATE_GENES(
            POLCA.out.polished,
            ch_ref_gff ?: [],
            ch_ref_fasta ?: [],
            ch_ref_gb ?: [],
            ch_trnascan ?: [],
            skip_trna,
        )
    }

    emit:
    polished   = POLCA.out.polished
    qc         = QC_SHORT.out.report
    asm_logs   = GETORGANELLE_ASM.out.logs
    vcf        = POLCA.out.vcf
    annotation = ANNOTATE_GENES.out
}
