//
// STRUCTURAL VARIATION
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GRIDSS as GRIDSS } from '../../../modules/nf-core/gridss/gridss/main'

workflow STRUCTURAL_VARIATION {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index (tuple)
        ch_fasta     // channel: [fasta]
        ch_fasta_fai // channel: fai

    main:

    ch_versions = Channel.empty()

    // Only one of the following should be run
    GRIDSS(ch_reads, ch_map_index, sort)

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
