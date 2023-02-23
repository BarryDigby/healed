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

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAGE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

gtf                = params.gtf                ? Channel.fromPath(params.gtf).collect()                       : Channel.empty()
fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()           : Channel.empty()

// Fails when wrongfull extension for intervals file
if (params.wes) {
    if (params.intervals && !params.intervals.endsWith("bed")) exit 1, "Target file specified with `--intervals` must be in BED format for targeted data"
    else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
} else if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("list")) exit 1, "Intervals file must end with .bed, .list, or .interval_list"



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE GENOME STEPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Enforce that the tool names are valid prior to this i.e json schema, passing an incorrect one on its own returns 'Cannot get property 'key' on null object'

// set up key value pairs
def dna_aligners_hashmap = [
    'bwa':'abra2,mutect2'
]

// combine all possible DNA tools that require aligners
dna_snv_indel_list = params.dna_snv_indel ? params.dna_snv_indel.split(',').collect{ it.trim().toLowerCase() } : []
dna_cnv_list = params.dna_cnv ? params.dna_cnv.split(',').collect{ it.trim().toLowerCase() } : []

dna_tools = dna_snv_indel_list + dna_cnv_list

def prepareDNAIndex = []
if(dna_tools.size() > 0) {
    for ( tool in dna_tools ) {
        prepareDNAIndex << dna_aligners_hashmap.find{ it.value.contains( tool ) }.key
    }
}

prepareDNAIndex.unique { a, b -> a <=> b }

def rna_aligners_hashmap = [
    'star':'star_salmon',
    'star_fusion':'star_fusion',
    'arriba':'arriba'
]

// combine all possible DNA tools that require aligners
rna_quant_list = params.rna_quant ? params.rna_quant.split(',').collect{ it.trim().toLowerCase() } : []
rna_fusion_list = params.rna_fusion ? params.rna_fusion.split(',').collect{ it.trim().toLowerCase() } : []

rna_tools = rna_quant_list + rna_fusion_list

def prepareRNAIndex = []
if(rna_tools.size() > 0) {
    for ( tool in rna_tools ) {
        prepareRNAIndex << rna_aligners_hashmap.find{ it.value.contains( tool ) }.key
    }
}

prepareRNAIndex.unique { a, b -> a <=> b }

// combine both lists to trigger genome index generation

prepareGenomeIndex = prepareRNAIndex + prepareDNAIndex

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


// Use subworkflows for everything!
include { BAM_ABRA2                   } from '../subworkflows/local/bam_abra2'
include { BAM_MARKDUPLICATES          } from '../subworkflows/local/bam_markduplicates'
include { FASTQ_FASTQC_FASTP          } from '../subworkflows/local/fastq_fastqc_fastp'
include { FASTQ_ALIGN_DNA             } from '../subworkflows/local/fastq_align_dna'
include { FASTQ_ALIGN_RNA             } from '../subworkflows/local/fastq_align_rna'
include { PREPARE_GENOME              } from '../subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS           } from '../subworkflows/local/prepare_intervals'
include { QUANTIFY_SALMON             } from '../subworkflows/local/quantify_salmon'
include { STARFUSION                  } from '../subworkflows/local/starfusion'
//include { STRUCTURAL_VARIATION        } from '../subworkflows/local/structural_variation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def multiqc_report = []

workflow HEALED {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                PARSE INPUT SAMPLESHEET
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Parse the input samplesheet.csv file to collect metadata, converting to channel con-
    taining metadata and reads. Sample output given below.

    Subworkflow, module files:
    - subworkflows/local/input_check
        - modules/local/samplesheet_check

    Parameters                                                               Explanation
    - params.fasta                                                             RG header
    - params.seq_center                                                        RG header
    - params.seq_platform                                                      RG header
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    INPUT_CHECK(
        ch_input
    )

    // Gather versions
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT_CHECK.out.reads.view()
    [[patient:HCC1395, assay:dna, status:normal, sample:HCC1395N_T1, lane:1, id:HCC1395N_T1-1, numLanes:1, read_group:"@RG\tID:null.HCC1395N_T1.1\tPU:1\tSM:HCC1395_HCC1395N_T1\tLB:HCC1395N_T1\tDS:null\tPL:ILLUMINA", data_type:fastq, size:1], [/data/github/healed/dev_data/HCC1395N_R1.fastq.gz, /data/github/healed/dev_data/HCC1395N_R2.fastq.gz]]
    [[patient:HCC1395, assay:dna, status:tumor,  sample:HCC1395T_T1, lane:1, id:HCC1395T_T1-1, numLanes:1, read_group:"@RG\tID:null.HCC1395T_T1.1\tPU:1\tSM:HCC1395_HCC1395T_T1\tLB:HCC1395T_T1\tDS:null\tPL:ILLUMINA", data_type:fastq, size:1], [/data/github/healed/dev_data/HCC1395T_R1.fastq.gz, /data/github/healed/dev_data/HCC1395T_R2.fastq.gz]]
    [[patient:GM12878, assay:rna, status:normal, sample:GM12878N_T1, lane:1, id:GM12878N_T1-1, numLanes:1, read_group:"@RG\tID:null.GM12878N_T1.1\tPU:1\tSM:GM12878_GM12878N_T1\tLB:GM12878N_T1\tDS:null\tPL:ILLUMINA", data_type:fastq, size:1, strandedness:reverse], [/data/github/healed/dev_data/GM12878_R1.fastq.gz, /data/github/healed/dev_data/GM12878_R2.fastq.gz]]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                FASTQC, FASTP, SPLIT_FASTQ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Perform FASTQC on raw reads, perform read trimming and split fastq files to paralle-
    lize downstream alignment steps.

    Subworkflow, module files:
    - subworkflows/local/fastq_fastqc_fastp
        - modules/nf-core/fastqc
        - modules/nf-core/fastp

    Config file:
    - conf/modules/fastq_fastqc_fastp.config

    Parameters                                                               Explanation
    - params.skip_fastqc                                                    Skip FASTQC?
    - params.trim_fastq                                          Trim FASTQ using FASTP?
    - params.split_fastq                                        Split FASTQ using FASTP?
    - params.save_trimmed                                        Save FASTP output file?
    - params.save_split_fastqs                                  Save split FASTP FASTQs?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_FASTQC_FASTP(
        INPUT_CHECK.out.reads
    )

    // Gather versions
    ch_versions = ch_versions.mix(FASTQ_FASTQC_FASTP.out.versions)

    // Gather reports
    ch_reports  = ch_reports.mix(FASTQ_FASTQC_FASTP.out.reports)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ_FASTQC_FASTP.out.reads.view()
    params.split_fastq 500000
    [[assay:dna, data_type:fastq, id:HCC1395T_T1-1, numLanes:1, patient:HCC1395, read_group:"@RG\tID:null.HCC1395T_T1.1\tPU:1\tSM:HCC1395_HCC1395T_T1\tLB:HCC1395T_T1\tDS:null\tPL:ILLUMINA", sample:HCC1395T_T1, size:2, status:tumor], [0001.HCC1395T_T1-1_1.fastp.fastq.gz, 0001.HCC1395T_T1-1_2.fastp.fastq.gz]]
    [[assay:dna, data_type:fastq, id:HCC1395T_T1-1, numLanes:1, patient:HCC1395, read_group:"@RG\tID:null.HCC1395T_T1.1\tPU:1\tSM:HCC1395_HCC1395T_T1\tLB:HCC1395T_T1\tDS:null\tPL:ILLUMINA", sample:HCC1395T_T1, size:2, status:tumor], [0002.HCC1395T_T1-1_1.fastp.fastq.gz, 0002.HCC1395T_T1-1_2.fastp.fastq.gz]]

    params.split_fastq 0
    [[patient:HCC1395, assay:dna, status:tumor, sample:HCC1395T_T1, lane:1, id:HCC1395T_T1-1, numLanes:1, read_group:"@RG\tID:null.HCC1395T_T1.1\tPU:1\tSM:HCC1395_HCC1395T_T1\tLB:HCC1395T_T1\tDS:null\tPL:ILLUMINA", data_type:fastq, size:1], [HCC1395T_T1-1_1.fastp.fastq.gz, HCC1395T_T1-1_2.fastp.fastq.gz]]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                PREPARE INTERVALS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    For DNA post alignment steps (i.e GATK best practices or Variant Calling) it is use-
    ful to perform scatter-gather on BAM files to parellelize processes. For WES, it is
    expected the user provides the capture kit BED file. For WGS, the fasta_fai file is
    used to derive a per-centromere intervals file.

    Subworkflow, module files:
    - subworkflows/local/prepare_intervals
        - modules/local/build_intervals
        - modules/local/create_intervals_bed
        - modules/nf-core/gatk4/intervallisttobed
        - modules/nf-core/samtools/faidx
        - modules/nf-core/tabix/bgziptabix

    Config file:
    - conf/modules/prepare_intervals.config

    Parameters                                                               Explanation
    - params.intervals                                   Intervals file provided by user
    - params.nucleotides_per_second    Estimate chunk size time: end-start / nuc_per_sec
    - params.no_intervals                 Output Channel.value([]) in place of intervals
    - params.save_reference                                Save generated interval files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    PREPARE_INTERVALS(
        fasta
    )

    // Gather versions
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    // Collect created files
    fasta_fai                   = PREPARE_INTERVALS.out.fasta_fai
    intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_for_preprocessing = params.wes          ? intervals_bed_combined : Channel.value([])
    intervals                   = PREPARE_INTERVALS.out.intervals_bed
    intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Tested with full Twist exome BED file to delineate channels
    WES, intervals file provided:
    chr21.fa.fai
    [twist.bed]
    [twist.bed]
    [chr5_141344034-141346523.bed, 61] ...
    [[chr5_141344034-141346523.bed.gz, chr5_141344034-141346523.bed.gz.tbi], 61] ...

    WGS, no interval file provided:
    chr21.fa.fai
    [chr21.fa.bed]
    []
    [chr21_1-46709983.bed, 1]
    [[chr21_1-46709983.bed.gz, chr21_1-46709983.bed.gz.tbi], 1]

    Explained:
    Fasta fai file
    All intervals in one file, or empty value channel
    All intervals for proprocessing. WGS can be empty, MOSDEPTH does not need them
    Sharded interval files with total number intervals
    Gzipped sharded interval files + index with total number intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                PREPARE GENOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Using Groovy hashmaps (think python dictionaries) parse the analysis tools selected
    by the user. From here we can work backwards and assign values to a list that defi-
    ne which genome index files must be prepared for downstream alignment. (See PREPARE
    GENOME STEPS in preamble).

    Subworkflow, module files:
    - subworkflows/local/prepare_genome
        - modules/nf-core/bwa/index
        - modules/nf-core/star/genomegenerate

    Config file:
    - conf/modules/prepare_genome.config

    Parameters                                                               Explanation
    - params.arriba_blacklist                                   Path to Arriba blacklist
    - params.arriba_known_fusions                           Path to Arriba known fusions
    - params.arriba_protein_domains                       Path to Arriba protein domains
    - params.bwa_index                                       Path to BWA index directory
    - params.star_index                                    Path to STAR inbdex directory
    - params.no_intervals                 Output Channel.value([]) in place of intervals
    - params.save_reference                                Save generated interval files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    PREPARE_GENOME(
        fasta,
        fasta_fai,
        gtf,
        prepareGenomeIndex
    )

    // Gather versions
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_GENOME.out.bwa.view()
    [[id:fasta.baseName], bwa]

    PREPARE_GENOME.out.star.view()
    star/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                SPLIT ASSAYS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    After processing reads, it makes sense to split the channel based on genomic assays.
    To my knowledge, no tool takes as input both DNA and RNA sequencing files simulaten-
    eously.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_reads = FASTQ_FASTQC_FASTP.out.reads
    ch_reads.branch{ meta, reads ->
        dna: meta.assay == 'dna'; return [meta, reads]
        rna: meta.assay == 'rna'; return [meta, reads]
    }.set{ ch_reads_to_map }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                FASTQ ALIGN DNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Channel processing is performed inside the workflow whereby the meta.id key is upda-
    ted to meta.sample if the sample is split across multiple lanes or fastq files such
    that they can be grouped after alignment using bam_merge_index subworkflow. A lot of
    the nuance with this subworkflow is based around how the user wants to save output
    files i.e as BAM or CRAM. The previously generated list 'dna_tools' can be used to
    trigger process events - much like prepare genome.

    Subworkflow, module files:
    - subworkflows/local/fastq_align_dna
        - modules/nf-core/bwa/mem
        - subworkflows/local/bam_merge_index
            - modules/nf-core/samtools/merge
            - modules/nf-core/samtools/index

    Config file:
    - conf/modules/fastq_align_dna.config

    Parameters                                                               Explanation
    params.dna_snv_indel           Tools for SNV/INDEL detection. Informs aligner choice
    params.dna_sv        Tools for detecting Structural Varation. Informs aligner choice
    params.save_mapped                                                Save mapped files?
    params.save_output_as_bam                                  Save as BAM? If not, CRAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_ALIGN_DNA(
        ch_reads_to_map.dna,
        PREPARE_GENOME.out.bwa_index,
        fasta,
        fasta_fai,
        dna_tools
    )

    // Gather versions
    ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ_ALIGN_DNA.out.cram.view():
    [[id:HCC1395T_T1, data_type:bam, patient:HCC1395, sample:HCC1395T_T1, status:tumor], HCC1395T_T1.sorted.cram, HCC1395T_T1.sorted.cram.crai]
    [[id:HCC1395N_T1, data_type:bam, patient:HCC1395, sample:HCC1395N_T1, status:normal], HCC1395N_T1.sorted.cram, HCC1395N_T1.sorted.cram.crai]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                FASTQ ALIGN RNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Firstly, FASTQ files over multiple lanes must be consolidated prior to alignment. T-
    -his is the opposite behaviour to DNA alignemnt which can be merged post alignment.
    Once alignment has been completed, the BAM files are sorted (if necessary), indexed
    and subject to stats, idxstats, flagstats. Saved as BAM or CRAM.

    STAR 'Quantification' is for downstream use by Salmon. We are interested only in the
    toTranscriptome.out.bam file for this step.

    Subworkflow, module files:
    - subworkflows/local/fastq_align_rna
        - modules/nf-core/cat/fastq
        - modules/nf-core/star/align
        - subworkflows/nf-core/bam_sort_stats_samtools
        - modules/nf-core/samtools/convert

    Config file:
    - conf/modules/fastq_align_rna.config

    Parameters                                                               Explanation
    params.rna_quant            Tool used for RNA quantification. Informs aligner choice
    params.rna_fusion         Tool used for RNA-Fusion detection. Informs aligner choice
    params.save_mapped                                                Save mapped files?
    params.save_output_as_bam                                  Save as BAM? If not, CRAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    FASTQ_ALIGN_RNA(
        ch_reads_to_map.rna,
        PREPARE_GENOME.out.star_index,
        fasta,
        fasta_fai,
        gtf,
        rna_tools
    )

    // Gather verions,reports
    ch_versions = ch_versions.mix(FASTQ_ALIGN_RNA.out.versions)
    ch_reports  = ch_reports.mix(FASTQ_ALIGN_RNA.out.reports)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ_ALIGN_RNA.out.quant_bam_transcript.view():
    [[id:GM12878N, sample:GM12878N, patient:GM12878, status:normal, strandedness:reverse], GM12878N.quant.bam, GM12878N.quant.bam.bai]
    FASTQ_ALIGN_RNA.out.fusion_junctions.view():
    [[id:fusion2, sample:fusion2, patient:fusion2, status:tumor, strandedness:unstranded], fusion2.fusion.Chimeric.out.junction]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                SALMON QUANTIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The transcriptome bam generated by FASTQ_ALIGN_RNA.out.quant_bam_transcript is prov-
    ded to salmon to peform quantification (based on alignments). There is some polishi-
    ng performed on the input GTF file whereby you can specify the attribute type to gr-
    oup features in the GTF file when running Salmon (i.e )

    Subworkflow, module files:
    - subworkflows/local/quantify_salmon
        - modules/nf-core/salmon/quant
        - modules/local/salmon_tx2gene
        - modules/local/salmon_tximport
        - modules/local/salmon_summarizedexperiment

    Config file:
    - conf/modules/quantify_salmon.config

    Parameters                                                               Explanation
    params.gtf_group_features                                  Second column for tx2gene
    params.gtf_extra_attributes                                 Third column for tx2gene
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    QUANTIFY_SALMON(
        FASTQ_ALIGN_RNA.out.salmon_bam_transcript.map{ meta, bam, bai -> [ meta, bam ] },
        ch_dummy_file,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.filter_gtf,
        true,
        params.salmon_quant_libtype ?: ''
    )

    // Gather versions
    ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Outputs are not yet consumed downstream, omitting channel structure until required.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                RNA FUSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utilise outputs from FASTQ_ALIGN_RNA specific to RNA Fusion. STAR-FUSION requires as
    input the junctions file output by STAR. Bring the sequencing reads too for fusion-
    inspector?

    Subworkflow, module files:
    - subworkflows/local/star_fusion
        - modules/local/starfusion/detect

    Config file:
    - conf/modules/starfusion.config

    Parameters                                                               Explanation

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    STARFUSION(
        FASTQ_ALIGN_RNA.out.starfusion_junctions,
        PREPARE_GENOME.out.ctat_genome_lib,
        rna_tools
    )

    // Gather versions
    ch_versions = ch_versions.mix(STARFUSION.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_MARKDUPLICATES.out.cram.view()
    [[id:HCC1395N_T1, data_type:bam, patient:HCC1395, sample:HCC1395N_T1, status:normal], HCC1395N_T1.md.cram, HCC1395N_T1.md.cram.crai]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                MARK DUPLICATES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Use aligned CRAMs as input, perform GATK MarkDuplicates in addition to samtools sta-
    ts and mosdepth stats. In and out, pretty straight forward.

    Subworkflow, module files:
    - subworkflows/local/bam_markduplicates
        - modules/nf-core/gatk4/markduplicates
        - modules/nf-core/samtools/index
        - subworkflows/local/cram_qc_mosdepth_samtools
            - modules/nf-core/samtools/stats
            - modules/nf-core/mosdepth

    Config file:
    - conf/modules/markduplicates.config

    Parameters                                                               Explanation
    params.rna_quant            Tool used for RNA quantification. Informs aligner choice
    params.rna_fusion         Tool used for RNA-Fusion detection. Informs aligner choice
    params.save_mapped                                                Save mapped files?
    params.save_output_as_bam                                  Save as BAM? If not, CRAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    BAM_MARKDUPLICATES(
        FASTQ_ALIGN_DNA.out.cram,
        fasta,
        fasta_fai,
        intervals_bed_combined
    )

    // Gather versions
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_MARKDUPLICATES.out.cram.view()
    [[id:HCC1395N_T1, data_type:bam, patient:HCC1395, sample:HCC1395N_T1, status:normal], HCC1395N_T1.md.cram, HCC1395N_T1.md.cram.crai]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*

    //
    // SAVE BAM AS CRAM
    //

    if (params.save_mapped) {

        // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
        BAM_MERGE_INDEX_SAMTOOLS(ch_bam_mapped)
        BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, ch_fai)

        // Gather used softwares versions
        ch_versions = ch_versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
        ch_versions = ch_versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
    }

    //
    // MARK DUPLICATES
    //

    // Use output from ALIGN_DNA [channel] tuple val(), path(bam)
    ch_for_markduplicates = ch_bam_mapped


    BAM_MARKDUPLICATES(
        ch_for_markduplicates,
        fasta,
        fasta_fai,
        intervals_for_preprocessing
    )

    // Capture output CRAM,CRAI file from MarkDuplicates
    ch_cram_markduplicates = BAM_MARKDUPLICATES.out.cram // [channel] cram, crai

    // Gather QC reports
    ch_reports  = ch_reports.mix(BAM_MARKDUPLICATES.out.qc.collect{meta, report -> report})

    // Gather used softwares versions
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES.out.versions)

    //
    // ABRA2 REALIGNMENT
    //

    if(params.abra2_realignment || dna_snv_indel_list.contains('abra2')) {

        BAM_ABRA2(
            ch_cram_markduplicates,
            fasta,
            ch_fai,
            gtf,
            junctions,
            intervals_for_preprocessing
        )

        // Gather versions
        ch_versions = ch_versions.mix(BAM_ABRA2.out.versions)

        // Gather ABRA2 crams
        ch_cram_abra2 = BAM_ABRA2.out.cram

    } else {
        ch_cram_abra2 = ch_cram_markduplicates
    }

    // Naming convention
    ch_for_recalibrator = ch_cram_abra2

    // take care of bams i.e merging if step skipped (hint place it in front of ABRA2 using mix on input chan , declare output from if statwement)
    // take care of bam to cram tracing too baz
    // make better comments on monday as you can't come back to this kind of shite after a break.

    //
    // STRUCTURAL VARIANTS
    //

   ///STRUCTURAL_VARIATION(
    //    ch_abrA
    //) */

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
    //ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))

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
