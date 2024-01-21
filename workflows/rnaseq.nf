/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['bowtie2'],
    trimmers       : ['fastp'],
    pseudoaligners : ['salmon'],
    // rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnaseq.initialise(workflow, params, log, valid_params)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.star_index, params.salmon_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc/main'
include { DESEQ2_QC as DESEQ2_QC_SALMON      } from '../modules/local/deseq2_qc/main'
// include { DUPRADAR                           } from '../modules/local/dupradar'
include { MULTIQC                            } from '../modules/local/multiqc/main'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_BOWTIE2     } from '../subworkflows/local/align_bowtie2'
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
// include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon/main'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow RNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.bowtie2_index,
        params.salmon_index,
        params.gencode,
        biotype,
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
            [ meta + [id: new_id], fastq ]
    }
    .groupTuple()
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_cat_fastq
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    PREPARE_GENOME.out.fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.salmon_index,
        !params.salmon_index && !('salmon' in prepareToolIndices)
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: WorkflowRnaseq.getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

        FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
            .map { meta, zip_file -> [ meta ] + WorkflowRnaseq.getTotalSequencesFastqc(workflow, meta, zip_file) }
            .set { ch_acc_reads_json }
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_trim_read_count
        .map {
            meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    ch_sortmerna_multiqc = Channel.empty()

    //
    // SUBWORKFLOW: Alignment with BOWTIE2 and gene/transcript quantification with Salmon
    //
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_bowtie2_multiqc            = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bowtie2_salmon') {
        ALIGN_BOWTIE2 (
            ch_filtered_reads,
            PREPARE_GENOME.out.bowtie2_index,
            PREPARE_GENOME.out.gtf,
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = ALIGN_BOWTIE2.out.bam
        ch_genome_bam_index  = ALIGN_BOWTIE2.out.bai
        ch_transcriptome_bam = ALIGN_BOWTIE2.out.bam_transcript
        ch_samtools_stats    = ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2.out.idxstats
        ch_star_multiqc      = ALIGN_BOWTIE2.out.log_final
        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_BOWTIE2.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions)

        // //
        // // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        // //
        // if (params.with_umi) {
        //     // Deduplicate genome BAM file before downstream analysis
        //     BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
        //         ch_genome_bam.join(ch_genome_bam_index, by: [0]),
        //         params.umitools_dedup_stats
        //     )
        //     ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
        //     ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
        //     ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
        //     ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
        //     ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
        //     if (params.bam_csi_index) {
        //         ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
        //     }
        //     ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)

        //     // Co-ordinate sort, index and run stats on transcriptome BAM
        //     BAM_SORT_STATS_SAMTOOLS (
        //         ch_transcriptome_bam,
        //         PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        //     )
        //     ch_transcriptome_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        //     ch_transcriptome_sorted_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

        //     // Deduplicate transcriptome BAM file before read counting with Salmon
        //     BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
        //         ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
        //         params.umitools_dedup_stats
        //     )

        //     // Name sort BAM before passing to Salmon
        //     SAMTOOLS_SORT (
        //         BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam
        //     )

        //     // Only run prepare_for_rsem.py on paired-end BAM files
        //     SAMTOOLS_SORT
        //         .out
        //         .bam
        //         .branch {
        //             meta, bam ->
        //                 single_end: meta.single_end
        //                     return [ meta, bam ]
        //                 paired_end: !meta.single_end
        //                     return [ meta, bam ]
        //         }
        //         .set { ch_umitools_dedup_bam }

        //     // Fix paired-end reads in name sorted BAM file
        //     // See: https://github.com/nf-core/rnaseq/issues/828
        //     UMITOOLS_PREPAREFORSALMON (
        //         ch_umitools_dedup_bam.paired_end
        //     )
        //     ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

        //     ch_umitools_dedup_bam
        //         .single_end
        //         .mix(UMITOOLS_PREPAREFORSALMON.out.bam)
        //         .set { ch_transcriptome_bam }
        // }

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        QUANTIFY_BOWTIE2_SALMON (
            ch_transcriptome_bam,
            ch_dummy_file,
            PREPARE_GENOME.out.transcript_fasta,
            PREPARE_GENOME.out.gtf,
            true,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(QUANTIFY_BOWTIE2_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_BOWTIE2_SALMON (
                QUANTIFY_BOWTIE2_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_BOWTIE2_SALMON.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_BOWTIE2_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_BOWTIE2_SALMON.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with BOWTIE2 and gene/transcript quantification with RSEM
    //
    ch_rsem_multiqc = Channel.empty()

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    ch_hisat2_multiqc = Channel.empty()

    //
    // Filter channels to get samples that passed BOWTIE2 minimum mapping percentage
    //
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_multiqc
            .map { meta, align_log -> [ meta ] + WorkflowRnaseq.getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bam_index
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam_index }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_mapped_reads[meta.id] = true
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    pass_mapped_reads[meta.id] = false
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
    }

    //
    // MODULE: Run Preseq
    //
    ch_preseq_multiqc = Channel.empty()

    //
        
    }

    // //
    // // MODULE: STRINGTIE
    // //
    // if (!params.skip_alignment && !params.skip_stringtie) {
    //     STRINGTIE_STRINGTIE (
    //         ch_genome_bam,
    //         PREPARE_GENOME.out.gtf
    //     )
    //     ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    // }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    ch_featurecounts_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {

        PREPARE_GENOME
            .out
            .gtf
            .map { WorkflowRnaseq.biotypeInGtf(it, biotype, log) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(PREPARE_GENOME.out.gtf)
            .combine(biotype_in_gtf)
            .filter { it[-1] }
            .map { it[0..<it.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_featurecounts_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    //
        )
    }

    //
    // MODULE: Downstream QC steps
    //
    ch_qualimap_multiqc           = Channel.empty()
    ch_dupradar_multiqc           = Channel.empty()
    ch_bamstat_multiqc            = Channel.empty()
    ch_inferexperiment_multiqc    = Channel.empty()
    ch_innerdistance_multiqc      = Channel.empty()
    ch_junctionannotation_multiqc = Channel.empty()
    ch_junctionsaturation_multiqc = Channel.empty()
    ch_readdistribution_multiqc   = Channel.empty()
    ch_readduplication_multiqc    = Channel.empty()
    ch_fail_strand_multiqc        = Channel.empty()
    ch_tin_multiqc                = Channel.empty()
    // if (!params.skip_alignment && !params.skip_qc) {
    //     if (!params.skip_qualimap) {
    //         QUALIMAP_RNASEQ (
    //             ch_genome_bam,
    //             PREPARE_GENOME.out.gtf
    //         )
    //         ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results
    //         ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
    //     }

    //     if (!params.skip_dupradar) {
    //         DUPRADAR (
    //             ch_genome_bam,
    //             PREPARE_GENOME.out.gtf
    //         )
    //         ch_dupradar_multiqc = DUPRADAR.out.multiqc
    //         ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
    //     }

    //     if (!params.skip_rseqc && rseqc_modules.size() > 0) {
    //         BAM_RSEQC (
    //             ch_genome_bam.join(ch_genome_bam_index, by: [0]),
    //             PREPARE_GENOME.out.gene_bed,
    //             rseqc_modules
    //         )
    //         ch_bamstat_multiqc            = BAM_RSEQC.out.bamstat_txt
    //         ch_inferexperiment_multiqc    = BAM_RSEQC.out.inferexperiment_txt
    //         ch_innerdistance_multiqc      = BAM_RSEQC.out.innerdistance_freq
    //         ch_junctionannotation_multiqc = BAM_RSEQC.out.junctionannotation_log
    //         ch_junctionsaturation_multiqc = BAM_RSEQC.out.junctionsaturation_rscript
    //         ch_readdistribution_multiqc   = BAM_RSEQC.out.readdistribution_txt
    //         ch_readduplication_multiqc    = BAM_RSEQC.out.readduplication_pos_xls
    //         ch_tin_multiqc                = BAM_RSEQC.out.tin_txt
    //         ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

    //         ch_inferexperiment_multiqc
    //             .map {
    //                 meta, strand_log ->
    //                     def inferred_strand = WorkflowRnaseq.getInferexperimentStrandedness(strand_log, 30)
    //                     pass_strand_check[meta.id] = true
    //                     if (meta.strandedness != inferred_strand[0]) {
    //                         pass_strand_check[meta.id] = false
    //                         return [ "$meta.id\t$meta.strandedness\t${inferred_strand.join('\t')}" ]
    //                     }
    //             }
    //             .collect()
    //             .map {
    //                 tsv_data ->
    //                     def header = [
    //                         "Sample",
    //                         "Provided strandedness",
    //                         "Inferred strandedness",
    //                         "Sense (%)",
    //                         "Antisense (%)",
    //                         "Undetermined (%)"
    //                     ]
    //                     WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
    //             }
    //             .set { ch_fail_strand_multiqc }
    //     }
    // }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_salmon_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon') {
        QUANTIFY_SALMON (
            ch_filtered_reads,
            PREPARE_GENOME.out.salmon_index,
            ch_dummy_file,
            PREPARE_GENOME.out.gtf,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_salmon_multiqc = QUANTIFY_SALMON.out.results
        ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_SALMON (
                QUANTIFY_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_SALMON.out.pca_multiqc
            ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_SALMON.out.versions)
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_acc_reads_json
            .mix(ch_acc_duplicates_json)
            .map { meta, json_file -> json_file }
            .collect()
            .map { files -> WorkflowRnaseq.formatAccJsonFiles(workflow, files) }
            .set { ch_acc_metrics }

        methods_description    = WorkflowRnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),
            ch_multiqc_logo.collect().ifEmpty([]),
            ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]),
            ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv').ifEmpty([]),
            ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv').ifEmpty([]),
            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),
            ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
            ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_aligner_pca_multiqc.collect().ifEmpty([]),
            ch_aligner_clustering_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
            ch_tin_multiqc.collect{it[1]}.ifEmpty([]),
            ch_acc_metrics.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}