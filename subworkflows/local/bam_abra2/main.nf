//
// ABRA2 REALIGNMENT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { ABRA2                                          } from '../../../modules/local/abra2/abra2/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_ABRA2          } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_ABRA2          } from '../../../modules/nf-core/samtools/convert/main'

workflow BAM_ABRA2 {
    take:
        cram                          // channel: [mandatory] meta, cram, crai
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        gtf                           // channel: [mandatory] gtf
        junctions                     // channel: [optional]  SJ.out.tab
        intervals                     // channel: [mandatory] intervals_preprocessing

    main:

    ch_versions = Channel.empty()

    // Input is CRAM, convert to BAM
        CRAM_TO_BAM_ABRA2(
            cram,
            fasta,
            fasta_fai
        )

    // Gather versions
    ch_versions = ch_versions.mix(CRAM_TO_BAM_ABRA2.out.versions)

    // Split normal/tumor BAM for ABRA2
    CRAM_TO_BAM_ABRA2.out.alignment_index.branch{ meta, bam, bai ->
                                            tumor:  meta.status == 'tumor'; return [meta, bam, bai]
                                            normal: meta.status =='normal'; return [meta, bam, bai]
                                        }.set{ ch_abra2_input }

    //empty tuple for abra2 = Channel.of([[], [], []]) // this works but you lose meta info.
    // set ifEmpty rules.
    abra2_normal = ch_abra2_input.normal.ifEmpty(Channel.of([ [], [], [] ]))
    abra2_tumor  = ch_abra2_input.tumor.ifEmpty(Channel.of([ [], [], [] ]))

    ABRA2(
        abra2_normal,
        abra2_tumor,
        fasta,
        gtf,
        [], // STAR junctions or []
        intervals,
        false
    )

    // Gather versions
    ch_versions = ch_versions.mix(ABRA2.out.versions)

    // Mix normal/tumor samples for downstream (one per line)
    ch_abra2_bams = ABRA2.out.abra2_normal.mix(ABRA2.out.abra2_tumor)

    // Convert ABRA2 BAM to CRAM
    BAM_TO_CRAM_ABRA2(
        ch_abra2_bams,
        fasta,
        fasta_fai
    )

    // Gather versions
    ch_versions = ch_versions.mix(BAM_TO_CRAM_ABRA2.out.versions)

    cram_abra2 = BAM_TO_CRAM_ABRA2.out.alignment_index

    emit:
        cram     = cram_abra2
        versions = ch_versions // channel: [ versions.yml ]
}
