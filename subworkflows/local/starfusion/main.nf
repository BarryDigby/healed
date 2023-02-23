//
// STARFUSION
//

include { STARFUSION as DETECT_FUSIONS       } from '../../../modules/local/starfusion/detect/main'

workflow STARFUSION {
    take:
        junctions
        ctat_genome_lib
        rna_tools

    main:

    ch_versions = Channel.empty()
    ch_reports = Channel.empty()

    //
    // STARFUSION
    //
    ch_fusion_predictions   = Channel.empty()
    ch_fusion_abridged      = Channel.empty()
    ch_fusion_coding_effect = Channel.empty()
    if('star_fusion' in rna_tools) {
        DETECT_FUSIONS (
            junctions,
            ctat_genome_lib
        )
        ch_fusion_predictions   = DETECT_FUSIONS.out.fusions
        ch_fusion_abridged      = DETECT_FUSIONS.out.abridged
        ch_fusion_coding_effect = DETECT_FUSIONS.out.coding_effect
        ch_versions             = ch_versions.mix(DETECT_FUSIONS.out.versions)
    }

    emit:
        fusion_predictions   = ch_fusion_predictions
        fusion_abridged      = ch_fusion_abridged
        fusion_coding_effect = ch_fusion_coding_effect
        versions             = ch_versions
        reports              = ch_reports
}
