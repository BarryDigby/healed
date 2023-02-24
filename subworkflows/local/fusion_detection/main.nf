//
// FUSION DETECTION
//

include { ARRIBA     } from '../../../modules/nf-core/arriba/main'
include { STARFUSION } from '../../../modules/local/starfusion/detect/main'

workflow FUSION_DETECTION {
    take:
        starfusion_junctions
        ctat_genome_lib
        fasta
        gtf
        arriba_bam
        blacklist
        known_fusions
        protein_domains
        rna_tools

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // ARRIBA
    //
    ch_arriba_fusions           = Channel.empty()
    ch_arriba_fusions_discarded = Channel.empty()
    if('arriba' in rna_tools) {
        ARRIBA(
            arriba_bam.map{meta, bam, bai -> [meta, bam]},
            fasta,
            gtf,
            blacklist,
            known_fusions,
            [], // structural variants, GRIDSS?
            [], // tags
            protein_domains
        )
        ch_arriba_fusions           = ARRIBA.out.fusions
        ch_arriba_fusions_discarded = ARRIBA.out.fusions_fail
        ch_versions                 = ch_versions.mix(ARRIBA.out.versions)
    }

    //
    // STARFUSION
    //
    ch_starfusion_predictions   = Channel.empty()
    ch_starfusion_abridged      = Channel.empty()
    ch_starfusion_coding_effect = Channel.empty()
    if('star_fusion' in rna_tools) {
        STARFUSION (
            starfusion_junctions,
            ctat_genome_lib
        )
        ch_starfusion_predictions   = DETECT_FUSIONS.out.fusions
        ch_starfusion_abridged      = DETECT_FUSIONS.out.abridged
        ch_starfusion_coding_effect = DETECT_FUSIONS.out.coding_effect
        ch_versions                 = ch_versions.mix(DETECT_FUSIONS.out.versions)
    }

    emit:
        arriba_fusions           = ch_arriba_fusions
        arriba_fusions_discarded = ch_arriba_fusions_discarded
        starfusion_predictions   = ch_starfusion_predictions
        starfusion_abridged      = ch_starfusion_abridged
        starfusion_coding_effect = ch_starfusion_coding_effect
        versions                 = ch_versions
        reports                  = ch_reports
}
