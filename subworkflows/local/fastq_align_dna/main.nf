//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'


workflow FASTQ_ALIGN_DNA {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index
        sort         // boolean: [mandatory] true -> sort, false -> don't sort

    main:

    ch_versions = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(ch_reads, ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort)

    // Get the bam files from the aligner
    // Only one aligner is run
    ch_bam_mapped = Channel.empty()
    ch_bam_mapped = ch_bam_mapped.mix(BWAMEM1_MEM.out.bam)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions.first())

    emit:
        bam      = ch_bam_mapped // channel: [ [meta], bam ]
        versions = ch_versions   // channel: [ versions.yml ]
}
