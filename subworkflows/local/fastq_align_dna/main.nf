//
// MAPPING
//

include { BWA_MEM                                 } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING } from '../../../modules/nf-core/samtools/convert/main'
include { BAM_MERGE_INDEX_SAMTOOLS                } from '../bam_merge_index_samtools'


workflow FASTQ_ALIGN_DNA {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index (tuple)
        fasta
        fasta_fai
        dna_tools

    main:

    ch_versions = Channel.empty()

    // If one lane per sample or one (set of) FASTQ
    // files per sample, change meta.id to meta.sample
    // as there is no need to merge the sample downstream
    ch_dna_reads_to_map = ch_reads.map{ meta, reads ->
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

    ch_bam_mapped = Channel.empty()
    sort = true // sort bam
    BWA_MEM(
        ch_dna_reads_to_map,
        ch_map_index,
        sort
    )
    ch_bam_mapped = ch_bam_mapped.mix(BWA_MEM.out.bam)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // Grouping the bams from the same samples not to stall the workflow
    ch_bams_to_merge = ch_bam_mapped.map{ meta, bam ->
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

    // Merge and Index subworkflow:
    // 1. Merge only occurs if meta.size is > 1
    // 2. Indexing occurs regardless
    BAM_MERGE_INDEX_SAMTOOLS(
        ch_bams_to_merge
    )

    // GATK preproc tools accept CRAM, so let's stick with
    // this file format going forward where applicable.
    // Recall the only reason to switch to BAM is for file
    // saving, or if a tool requires it.
    BAM_TO_CRAM_MAPPING(
        BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai,
        fasta,
        fasta_fai
    )

    // Gather used softwares versions
    ch_versions = ch_versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
    ch_versions = ch_versions.mix(BAM_TO_CRAM_MAPPING.out.versions)

    emit:
        cram     = BAM_TO_CRAM_MAPPING.out.alignment_index // channel: [ [meta], cram, crai ]
        versions = ch_versions   // channel: [ versions.yml ]
}
