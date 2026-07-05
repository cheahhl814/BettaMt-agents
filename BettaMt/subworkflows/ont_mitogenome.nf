// ===========================================================================
// ONT_MITOGENOME — subworkflow for Oxford Nanopore long-read assembly
//
//   reads, ref_mito          ─► BAIT_MITO        ─► flye_in
//   flye_in, size            ─► FLYE             ─► flye_gfa  ─►  GETORGANELLE_FILTER ─► clean
//   clean                    ─► FIRSTGENE        ─► trnF
//   clean, trnF              ─► CIRCLATOR        ─► rotated
//   rotated, baited_reads    ─► RACON            ─► polished (final mitogenome)
//
// All processes are imported from modules/local/*.nf and run via pixi.
// ===========================================================================
nextflow.enable.dsl = 2

include { BAIT_MITO          } from '../modules/local/bait_mito.nf'
include { FLYE               } from '../modules/local/flye.nf'
include { GETORGANELLE_FILTER} from '../modules/local/getorganelle_filter.nf'
include { FIRSTGENE          } from '../modules/local/firstgene.nf'
include { CIRCLATOR          } from '../modules/local/circlator.nf'
include { RACON              } from '../modules/local/racon.nf'
include { ANNOTATE_GENES     } from '../modules/local/annotate_genes.nf'

workflow ONT_MITOGENOME {

    take:
    ch_reads        // channel: [ path(fastq) ]
    ch_ref_mito     // channel: [ path(fasta) ]
    ch_ref_gff      // channel: [ path(gff3) ]  — optional, for gene annotation
    ch_ref_fasta    // channel: [ path(fasta) ] — optional, companion to ref_gff
    ch_ref_gb       // channel: [ path(gb) ]    — optional, alternative to gff+fasta
    ch_trnascan     // channel: [ path(binary) ] — optional
    skip_trna       // value:    boolean
    size            // value:    estimated genome size (e.g. "16k")
    taxon           // value:    tRNAscan-SE model (e.g. "vertebrate")

    main:
    BAIT_MITO(ch_reads, ch_ref_mito)
    FLYE(BAIT_MITO.out.baited, size)
    GETORGANELLE_FILTER(FLYE.out.gfa)
    FIRSTGENE(GETORGANELLE_FILTER.out.contig, taxon)
    CIRCLATOR(GETORGANELLE_FILTER.out.contig, FIRSTGENE.out.fasta)
    RACON(CIRCLATOR.out.contig, BAIT_MITO.out.baited)

    // Gene annotation — only if a reference annotation was provided
    if (ch_ref_gff || ch_ref_gb) {
        ANNOTATE_GENES(
            RACON.out,
            ch_ref_gff ?: [],
            ch_ref_fasta ?: [],
            ch_ref_gb ?: [],
            ch_trnascan ?: [],
            skip_trna,
        )
    }

    emit:
    polished   = RACON.out
    circlator  = CIRCLATOR.out.contig
    baited     = BAIT_MITO.out.baited
    flye_gfa   = FLYE.out.gfa
    annotation = ANNOTATE_GENES.out
}
