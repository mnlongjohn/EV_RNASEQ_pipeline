//
// This file holds several functions specific to the workflow/rnaseq.nf in the nf-core/rnaseq pipeline
//

import nextflow.Nextflow
import java.util.zip.ZipFile
import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import groovy.text.SimpleTemplateEngine

class WorkflowRnaseq {

    //
    // Check and validate parameters
    //
    public static void initialise(workflow, params, log, valid_params) {
        if (!params.fasta) {
            Nextflow.error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
        }

        if (!params.gtf && !params.gff) {
            Nextflow.error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
        }

        if (params.gtf) {
            if (params.gff) {
                gtfGffWarn(log)
            }
            // if (params.genome == 'GRCh38' && params.gtf.contains('Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf')) {
            //     ncbiGenomeWarn(log)
            // }
            if (params.gtf.contains('/UCSC/') && params.gtf.contains('Annotation/Genes/genes.gtf')) {
                ucscGenomeWarn(log)
            }
        }

        if (params.transcript_fasta) {
            transcriptsFastaWarn(log)
        }

        if (!params.skip_trimming) {
            if (!valid_params['trimmers'].contains(params.trimmer)) {
                Nextflow.error("Invalid option: '${params.trimmer}'. Valid options for '--trimmer': ${valid_params['trimmers'].join(', ')}.")
            }
        }

        if (!params.skip_alignment) {
            if (!valid_params['aligners'].contains(params.aligner)) {
                Nextflow.error("Invalid option: '${params.aligner}'. Valid options for '--aligner': ${valid_params['aligners'].join(', ')}.")
            }
        } else {
            skipAlignmentWarn(log)
        }

        if (!params.skip_pseudo_alignment && params.pseudo_aligner) {
            if (!valid_params['pseudoaligners'].contains(params.pseudo_aligner)) {
                Nextflow.error("Invalid option: '${params.pseudo_aligner}'. Valid options for '--pseudo_aligner': ${valid_params['pseudoaligners'].join(', ')}.")
            } else {
                if (!(params.salmon_index || params.transcript_fasta || (params.fasta && (params.gtf || params.gff)))) {
                    Nextflow.error("To use `--pseudo_aligner 'salmon'`, you must provide either --salmon_index or --transcript_fasta or both --fasta and --gtf / --gff.")
                }
            }
        }

        // Warn if --additional_fasta provided with aligner index
        if (!params.skip_alignment && params.additional_fasta) {
            def index = ''
            if (params.aligner == 'star_salmon' && params.star_index) {
                index = 'star'
            }
            if (index) {
                additionaFastaIndexWarn(index, log)
            }
        }

        // // Check which RSeQC modules we are running
        // def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
        // if ((valid_params['rseqc_modules'] + rseqc_modules).unique().size() != valid_params['rseqc_modules'].size()) {
        //     Nextflow.error("Invalid option: ${params.rseqc_modules}. Valid options for '--rseqc_modules': ${valid_params['rseqc_modules'].join(', ')}")
        // }

        // initialise the custom acc json file
        def output_path = "${workflow.launchDir}/work/ACC_Guidelines"
        def output_dir = new File(output_path)
        // Check to insure that the dir exist before writing
        if (!output_dir.exists()) {
            output_dir.mkdirs()
        }
    }

    //
    // Function to check whether biotype field exists in GTF file
    //
    public static Boolean biotypeInGtf(gtf_file, biotype, log) {
        def hits = 0
        gtf_file.eachLine { line ->
            def attributes = line.split('\t')[-1].split()
            if (attributes.contains(biotype)) {
                hits += 1
            }
        }
        if (hits) {
            return true
        } else {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
                "  Biotype QC will be skipped to circumvent the issue below:\n" +
                "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
                "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            return false
        }
    }

    //
    // Function to generate an error if contigs in genome fasta file > 512 Mbp
    //
    public static void checkMaxContigSize(fai_file, log) {
        def max_size = 512000000
        fai_file.eachLine { line ->
            def lspl  = line.split('\t')
            def chrom = lspl[0]
            def size  = lspl[1]
            if (size.toInteger() > max_size) {
                def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
                    "  ${chrom}: ${size}\n\n" +
                    "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
                    "  Please see:\n" +
                    "  https://github.com/nf-core/rnaseq/issues/744\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.error(error_string)
            }
        }
    }

    //
    // Function that parses Salmon quant 'meta_info.json' output file to get inferred strandedness
    //
    public static String getSalmonInferredStrandedness(json_file) {
        def lib_type = new JsonSlurper().parseText(json_file.text).get('library_types')[0]
        def strandedness = 'reverse'
        if (lib_type) {
            if (lib_type in ['U', 'IU']) {
                strandedness = 'unstranded'
            } else if (lib_type in ['SF', 'ISF']) {
                strandedness = 'forward'
            } else if (lib_type in ['SR', 'ISR']) {
                strandedness = 'reverse'
            }
        }
        return strandedness
    }

    //
    // Function that parses TrimGalore log output file to get total number of reads after trimming
    //
    public static Integer getTrimGaloreReadsAfterFiltering(log_file) {
        def total_reads = 0
        def filtered_reads = 0
        log_file.eachLine { line ->
            def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
            def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
            if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
            if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
        return total_reads - filtered_reads
    }

    //
    // Function that parses fastqc's data output file to get the total number of reads
    //
    public static ArrayList getTotalSequencesFastqc(workflow, meta, zip_file) {
        def output_path = "${workflow.launchDir}/work/ACC_Guidelines"
        def output_file = new File("${output_path}/${meta.id}_fastqc_reads.json")
        
        // get total reads from zipped folder
        def total_reads = 0
        def zipFile = new ZipFile(zip_file[0].toString())
        zipFile.entries().each { entry ->
            if (!entry.isDirectory() && entry.name.contains('fastqc_data.txt')) {
                def inputStream = zipFile.getInputStream(entry) 
                def reader = new BufferedReader(new InputStreamReader(inputStream))
                reader.eachLine { line ->
                    def total_reads_matcher = line =~ /Total Sequences\s([\d\.]+)/
                    if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toInteger()
                }
            }
        }

        // print to json
        def data = [id: meta.id, total_reads: total_reads]
        output_file.withWriter { writer ->
            writer.write(JsonOutput.toJson(data))
        }
        // output_file.write(JsonOutput.toJson(data))
        def path_to_json = output_file.toPath()

        return [ path_to_json ]
    }

    //
    // Function that parses picard's markduplicate metrics output file to get percent duplication
    //
    public static ArrayList getPicardMarkDuplicatesRate(workflow, meta, log_file) {
        def output_path = "${workflow.launchDir}/work/ACC_Guidelines"
        def output_file = new File("${output_path}/${meta.id}_picard_duplicate.json")

        def percent_duplication = 0
        def line = log_file.readLines()[7]
        def metrics_list = line.split('\t')

        // if the metric exists print to json
        if (metrics_list) percent_duplication = metrics_list[8].toFloat() * 100
        def data = [id: meta.id, percent_duplication: percent_duplication]
        // output_file.write(JsonOutput.toJson(data))
        output_file.withWriter { writer ->
            writer.write(JsonOutput.toJson(data))
        }

        def path_to_json = output_file.toPath()

        return [ path_to_json ]
    }

    //
    // Function that parses and returns the alignment rate from the STAR log output
    //
    public static ArrayList getStarPercentMapped(params, align_log) {
        def percent_aligned = 0
        def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
        align_log.eachLine { line ->
            def matcher = line =~ pattern
            if (matcher) {
                percent_aligned = matcher[0][1].toFloat()
            }
        }

        def pass = false
        if (percent_aligned >= params.min_mapped_reads.toFloat()) {
            pass = true
        }
        return [ percent_aligned, pass ]
    }

    //
    // Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
    //
    public static ArrayList getInferexperimentStrandedness(inferexperiment_file, cutoff=30) {
        def sense        = 0
        def antisense    = 0
        def undetermined = 0
        inferexperiment_file.eachLine { line ->
            def undetermined_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
            def se_sense_matcher     = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
            def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
            def pe_sense_matcher     = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
            def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/
            if (undetermined_matcher) undetermined = undetermined_matcher[0][1].toFloat() * 100
            if (se_sense_matcher)     sense        = se_sense_matcher[0][1].toFloat() * 100
            if (se_antisense_matcher) antisense    = se_antisense_matcher[0][1].toFloat() * 100
            if (pe_sense_matcher)     sense        = pe_sense_matcher[0][1].toFloat() * 100
            if (pe_antisense_matcher) antisense    = pe_antisense_matcher[0][1].toFloat() * 100
        }
        def strandedness = 'unstranded'
        if (sense >= 100-cutoff) {
            strandedness = 'forward'
        } else if (antisense >= 100-cutoff) {
            strandedness = 'reverse'
        }
        return [ strandedness, sense, antisense, undetermined ]
    }

    //
    // Function merges the ACC json files together
    //
    public static ArrayList formatAccJsonFiles(workflow, file_array) {
        def output_path = "${workflow.launchDir}/work/ACC_Guidelines"
        def output_file = new File("${output_path}/ACC_RNASEQ_Guidelines_mqc.json")

        // read all the files and parse the json
        def jsonData = []

        file_array.each { file ->
            def data = new JsonSlurper().parseText(file.text)
            jsonData << data
        }

        // merge the data based on matching "id"
        def result = [:].withDefault { [] }
        jsonData.each { item ->
            result[item.id] << item
        }

        def mergedData = [:].withDefault { [] }
        result.each { id, items ->
            def mergedItem = [:]
            items.each { item ->
                mergedItem.putAll(item)
            }
            mergedData[id] = [total_reads: mergedItem.total_reads / 1000000, percent_duplication: mergedItem.percent_duplication]
        }

        // format table with data
        def outputJson = [
            id: "acc_rnaseq_guidelines",
            section_name: "ACC-MoH RNASEQ Metrics",
            description: "This table is a test to show if the samples run meet the ACC-MoH metrics of >80 million reads and <50% duplication.",
            plot_type: "table",
            pconfig: [
                id: "acc_rnaseq_guidelines_table",
                title: "ACC RNASEQ Guidelines Table",
            ],
            headers: [
                total_reads: [
                    title: "Total Reads",
                    description: "Total number of reads sequenced after filtering (millions)",
                    min: 0,
                    suffix: "M",
                    format: "{:,.0f}"
                ],
                percent_duplication: [
                    title: "Percent Duplication",
                    description: "Percentage of aligned reads that are duplicates",
                    max: 100,
                    min: 0,
                    suffix: "%",
                    format: "{:,.2f}"
                ]
            ],
            data: mergedData
        ]

        // write the json plot out
        output_file.withWriter { writer ->
            writer.write(JsonOutput.toJson(outputJson))
        }

        def path_to_json = output_file.toPath()

        return [ path_to_json ]
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Create MultiQC tsv custom content from a list of values
    //
    public static String multiqcTsvFromList(tsv_data, header) {
        def tsv_string = ""
        if (tsv_data.size() > 0) {
            tsv_string += "${header.join('\t')}\n"
            tsv_string += tsv_data.join('\n')
        }
        return tsv_string
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Print a warning if using GRCh38 assembly from igenomes.config
    //
    private static void ncbiGenomeWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  When using '--genome GRCh38' the assembly is from the NCBI and NOT Ensembl.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if using a UCSC assembly from igenomes.config
    //
    private static void ucscGenomeWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  When using UCSC assemblies the 'gene_biotype' field is absent from the GTF file.\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  If you would like to use the soft-masked Ensembl assembly instead please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/159#issuecomment-501184312\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if both GTF and GFF have been provided
    //
    private static void gtfGffWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Both '--gtf' and '--gff' parameters have been provided.\n" +
            "  Using GTF file as priority.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if using '--transcript_fasta'
    //
    private static void transcriptsFastaWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  '--transcript_fasta' parameter has been provided.\n" +
            "  Make sure transcript names in this file match those in the GFF/GTF file.\n\n" +
            "  Please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/753\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if --skip_alignment has been provided
    //
    private static void skipAlignmentWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  '--skip_alignment' parameter has been provided.\n" +
            "  Skipping alignment, genome-based quantification and all downstream QC processes.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning if using '--additional_fasta' and '--<ALIGNER>_index'
    //
    private static void additionaFastaIndexWarn(index, log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  When using '--additional_fasta <FASTA_FILE>' the aligner index will not\n" +
            "  be re-built with the transgenes incorporated by default since you have \n" +
            "  already provided an index via '--${index}_index <INDEX>'.\n\n" +
            "  Set '--additional_fasta <FASTA_FILE> --${index}_index false --gene_bed false --save_reference'\n" +
            "  to re-build the index with transgenes included and the index and gene BED file will be saved in\n" +
            "  'results/genome/index/${index}/' for re-use with '--${index}_index'.\n\n" +
            "  Ignore this warning if you know that the index already contains transgenes.\n\n" +
            "  Please see:\n" +
            "  https://github.com/nf-core/rnaseq/issues/556\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }
}
