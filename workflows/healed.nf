/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowHealed.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.bwa_index, params.star_index ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAGE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

gtf                = params.gtf                ? Channel.fromPath(params.gtf)                       : Channel.empty()
fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()                     : Channel.empty()

// Fails when wrongfull extension for intervals file
//if (params.wes) {
//    if (params.intervals && !params.intervals.endsWith("bed")) exit 1, "Target file specified with `--intervals` must be in BED format for targeted data"
//    else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
//} else if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("list")) exit 1, "Intervals file must end with .bed, .list, or .interval_list"



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE GENOME STEPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



// set up key value pairs
def dna_aligners_hashmap = [
    'bwa':'abra2,mutect2'
]

// combine all possible DNA tools that require aligners
dna_snv_indel_list = params.dna_snv_indel ? params.dna_snv_indel.split(',').collect{ it.trim().toLowerCase() } : []
dna_cnv_list = params.dna_cnv ? params.cnv.split(',').collect{ it.trim().toLowerCase() } : []

dna_tools = dna_snv_indel_list + dna_cnv_list

def prepareDNAIndex = []
if(dna_tools) {
    for ( tool in dna_tools ) {
        prepareDNAIndex << dna_aligners_hashmap.find{ it.value.contains( tool ) }.key
    }
}

prepareDNAIndex.unique { a, b -> a <=> b }

// same for rna, combine the two lists at the end.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                         } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM                } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING        } from '../modules/nf-core/samtools/convert/main'

include { BAM_MARKDUPLICATES          } from '../subworkflows/local/bam_markduplicates'
include { BAM_MERGE_INDEX_SAMTOOLS    } from '../subworkflows/local/bam_merge_index_samtools'
include { FASTQ_FASTQC_FASTP          } from '../subworkflows/local/fastq_fastqc_fastp'
include { FASTQ_ALIGN_DNA             } from '../subworkflows/local/fastq_align_dna'
include { PREPARE_GENOME              } from '../subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS           } from '../subworkflows/local/prepare_intervals'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def multiqc_report = []

workflow HEALED {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // Parse samplesheet
    //
    INPUT_CHECK(
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // FASTQC, FASTP, SPLIT FASTQ
    //
    FASTQ_FASTQC_FASTP(
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQ_FASTQC_FASTP.out.versions)
    ch_reports  = ch_reports.mix(FASTQ_FASTQC_FASTP.out.reports)

    //
    // PREPARE INTERVALS
    //
    PREPARE_INTERVALS(
        fasta
    )
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    // Collect created files
    fasta_fai                   = PREPARE_INTERVALS.out.fasta_fai
    intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined  // [interval.bed] all intervals in one file
    intervals_for_preprocessing = params.wes          ? intervals_bed_combined : []      // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals                   = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi

    //
    // PREPARE GENOME
    //
    PREPARE_GENOME(
        fasta,
        fasta_fai,
        gtf,
        prepareDNAIndex
    )

    //
    // SPLIT ASSAYS
    //

    ch_reads = FASTQ_FASTQC_FASTP.out.reads
    ch_reads.branch{ meta, reads ->
        dna: meta.assay == 'dna'; return [meta, reads]
        rna: meta.assay == 'rna'; return [meta, reads]
    }.set{ ch_reads_to_map }

    //
    // FASTQ ALIGN DNA
    //

    ch_dna_reads_to_map = ch_reads_to_map.dna.map{ meta, reads ->
        // update ID when no multiple lanes or splitted fastqs
        new_id = meta.size * meta.numLanes == 1 ? meta.sample : meta.id

        [[
            data_type:  meta.data_type,
            id:         new_id,
            numLanes:   meta.numLanes,
            patient:    meta.patient,
            read_group: meta.read_group,
            sample:     meta.sample,
            size:       meta.size,
            status:     meta.status,
            ],
        reads]
    }

    sort_bam = true
    FASTQ_ALIGN_DNA(
        ch_dna_reads_to_map,
        PREPARE_GENOME.out.bwa_index,
        sort_bam
    )

    //
    // ALIGNED DNA BAMS
    //

    // Grouping the bams from the same samples not to stall the workflow
    ch_bam_mapped = FASTQ_ALIGN_DNA.out.bam.map{ meta, bam ->
        numLanes = meta.numLanes ?: 1
        size     = meta.size     ?: 1

        // update ID to be based on the sample name
        // update data_type
        // remove no longer necessary fields:
        //   read_group: Now in the BAM header
        //     numLanes: Was only needed for mapping
        //         size: Was only needed for mapping
        new_meta = [
                    id:meta.sample,
                    data_type:"bam",
                    patient:meta.patient,
                    sample:meta.sample,
                    status:meta.status,
                ]

        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        [ groupKey(new_meta, numLanes * size), bam]
    }.groupTuple()

    //
    // SAVE BAM AS CRAM
    //

    if (params.save_mapped) {

        // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
        BAM_MERGE_INDEX_SAMTOOLS(ch_bam_mapped)
        BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)

        // Gather used softwares versions
        ch_versions = ch_versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
    }

    //
    // MARK DUPLICATES
    //

    ch_for_markduplicates = ch_bam_mapped
    intervals_for_preprocessing.view()
    BAM_MARKDUPLICATES(
        ch_for_markduplicates,
        fasta,
        fasta_fai,
        intervals_for_preprocessing
    )

    ch_cram_markduplicates = BAM_MARKDUPLICATES.out.cram

    // Gather QC reports
    ch_reports  = ch_reports.mix(BAM_MARKDUPLICATES.out.qc.collect{meta, report -> report})
    // Gather used softwares versions
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)

    //
    // COLLECT VERSIONS
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowHealed.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowHealed.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
