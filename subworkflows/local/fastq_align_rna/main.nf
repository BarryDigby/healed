//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { STAR_ALIGN as STAR_2PASS_BASIC } from '../modules/nf-core/star/align/main'

workflow FASTQ_ALIGN_RNA {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index (tuple)
        gtf          // channel: [mandatory] gtf

    main:

    ch_versions = Channel.empty()

    star_ignore_sjdbgtf = false // use genome GTF file
    seq_center          = params.seq_center ?: ''
    seq_platform        = 'ILLUMINA'

    STAR_2PASS_BASIC( ch_reads, ch_map_index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )

    ch_versions = ch_versions.mix(STAR_2PASS_BASIC.out.versions)

    // Get the bam files from the aligner
    ch_bam_mapped = Channel.empty()
    ch_junctions  = Channel.empty()
    ch_bam_mapped = ch_bam_mapped.mix(STAR_2PASS_BASIC.out.bam_sorted)
    ch_junctions  = ch_junctions.mix(STAR_2PASS_BASIC.out.junctions)

    emit:
        bam       = ch_bam_mapped // [channel] tuple
        junctions = ch_junctions  // [channel] tuple
        versions  = ch_versions   // [channel] path
}
