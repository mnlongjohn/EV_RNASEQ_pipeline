//
// Alignment with STAR
//

include { BOWTIE2_ALIGN          } from '../../modules/nf-core/bowtie2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow ALIGN_STAR_BOWTIE2 {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: /path/to/star/index/
    save_unaligned      // val
    sort_bam            // val
    fasta               // channel: /path/to/fasta

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with Bowtie2
    //
    ch_orig_bam       = Channel.empty()
    ch_log_final      = Channel.empty()
    ch_log_out        = Channel.empty()
    ch_log_progress   = Channel.empty()
    ch_bam_sorted     = Channel.empty()
    ch_bam_transcript = Channel.empty()
    ch_fastq          = Channel.empty()
    ch_tab            = Channel.empty()
    BOWTIE2_ALIGN ( reads, index, save_unaligned, sort_bam )
    ch_orig_bam       = BOWTIE2_ALIGN.out.bam
    ch_log_final      = BOWTIE2_ALIGN.out.log_final
    ch_log_out        = BOWTIE2_ALIGN.out.log_out
    ch_log_progress   = BOWTIE2_ALIGN.out.log_progress
    ch_bam_sorted     = BOWTIE2_ALIGN.out.bam_sorted
    ch_bam_transcript = BOWTIE2_ALIGN.out.bam_transcript
    ch_fastq          = BOWTIE2_ALIGN.out.fastq
    ch_tab            = BOWTIE2_ALIGN.out.tab
    ch_versions       = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( ch_orig_bam, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    orig_bam       = ch_orig_bam                    // channel: [ val(meta), bam            ]
    log_final      = ch_log_final                   // channel: [ val(meta), log_final      ]
    log_out        = ch_log_out                     // channel: [ val(meta), log_out        ]
    log_progress   = ch_log_progress                // channel: [ val(meta), log_progress   ]
    bam_sorted     = ch_bam_sorted                  // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript              // channel: [ val(meta), bam_transcript ]
    fastq          = ch_fastq                       // channel: [ val(meta), fastq          ]
    tab            = ch_tab                         // channel: [ val(meta), tab            ]

    bam            = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai            = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi            = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats          = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat       = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats       = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions       = ch_versions                    // channel: [ versions.yml ]
}
