//
// FUSION DETECTION
//

include { ARRIBA        } from '../../../modules/nf-core/arriba/main'
include { STARFUSION    } from '../../../modules/local/starfusion/detect/main'
include { FUSIONCATCHER } from '../../../modules/local/fusioncatcher/detect/main'

workflow FUSION_DETECTION {
    take:
        starfusion_junctions
        ctat_genome_lib
        fusioncatcher_lib
        reads
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
    // FUSIONCATCHER
    //
    ch_fusioncatcher_fusions = Channel.empty()
    ch_fusioncatcher_summary = Channel.empty()
    ch_fusioncatcher_log     = Channel.empty()
    if('fusioncatcher' in rna_tools) {
        FUSIONCATCHER(
            reads,
            fusioncatcher_lib
        )
        ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions
        ch_fusioncatcher_summary = FUSIONCATCHER.out.summary
        ch_fusioncatcher_log     = FUSIONCATCHER.out.log
        ch_versions              = ch_versions.mix(FUSIONCATCHER.out.versions)
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
        ch_starfusion_predictions   = STARFUSION.out.fusions
        ch_starfusion_abridged      = STARFUSION.out.abridged
        ch_starfusion_coding_effect = STARFUSION.out.coding_effect
        ch_versions                 = ch_versions.mix(STARFUSION.out.versions)
    }

    emit:
        arriba_fusions           = ch_arriba_fusions
        arriba_fusions_discarded = ch_arriba_fusions_discarded
        fusioncatcher_fusions    = ch_fusioncatcher_fusions
        fusioncatcher_summary    = ch_fusioncatcher_summary
        fusioncatcher_log        = ch_fusioncatcher_log
        starfusion_predictions   = ch_starfusion_predictions
        starfusion_abridged      = ch_starfusion_abridged
        starfusion_coding_effect = ch_starfusion_coding_effect
        versions                 = ch_versions
        reports                  = ch_reports
}
