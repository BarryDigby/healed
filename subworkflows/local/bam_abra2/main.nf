//
// ABRA2 REALIGNMENT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { ABRA2 as ABRA2_DNA                       } from '../../../modules/local/abra2/abra2/main'
include { ABRA2 as ABRA2_RNA                       } from '../../../modules/local/abra2/abra2/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM          } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM          } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_RNA      } from '../../../modules/nf-core/samtools/convert/main'

workflow BAM_ABRA2 {
    take:
        dna_cram                          // channel: meta, cram, crai
        rna_bam                           // channel: meta, bam, bai
        fasta                             // channel: fasta
        fasta_fai                         // channel: fasta_fai
        gtf                               // channel: gtf
        junctions                         // channel: SJ.out.tab
        intervals                         // channel: intervals_preprocessing
        dna_tools                         // [Array.list] for dev testing runs
        rna_tools                         // [Array.list] for dev testing runs

    main:

    // dev testing runs: when you want to test RNA or DNA and not both at once.

    ch_versions = Channel.empty()

    //
    // ABRA2 DNA
    //
    ch_abra2_dna = Channel.empty()
    if( dna_tools.size() > 0 && ( params.abra2_realignment || params.dna_snv_indel.contains('abra_cadabra') ) ){

        CRAM_TO_BAM(
            dna_cram,
            fasta,
            fasta_fai
        )
        ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions)

        CRAM_TO_BAM.out.alignment_index.branch{ meta, bam, bai ->
                        tumor:  meta.status == 'tumor'; return [meta, bam, bai]
                        normal: meta.status =='normal'; return [meta, bam, bai]
                    }.set{ ch_abra2_dna_input }
        abra2_dna_normal = ch_abra2_dna_input.normal.ifEmpty(Channel.of([ [], [], [] ]))
        abra2_dna_tumor  = ch_abra2_dna_input.tumor.ifEmpty(Channel.of([ [], [], [] ]))

        ABRA2_DNA(
            abra2_dna_normal,
            abra2_dna_tumor,
            fasta,
            [],
            [],
            intervals
        )

        ch_versions = ch_versions.mix(ABRA2_DNA.out.versions)

        BAM_TO_CRAM(
            ABRA2_DNA.out.abra2_normal.mix(ABRA2_DNA.out.abra2_tumor),
            fasta,
            fasta_fai
        )

        ch_abra2_dna = ch_abra2_dna.mix(BAM_TO_CRAM.out.alignment_index)
        ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)
    }

    //
    // ABRA2 RNA
    //
    ch_abra2_rna = Channel.empty()
    if( rna_tools.size() > 0 && params.abra2_realignment ){

        rna_bam.branch{ meta, bam, bai ->
                    tumor:  meta.status == 'tumor'; return [meta, bam, bai]
                    normal: meta.status =='normal'; return [meta, bam, bai]
        }.set{ ch_abra2_rna_input }
        abra2_rna_normal = ch_abra2_rna_input.normal.ifEmpty(Channel.of([ [], [], [] ]))
        abra2_rna_tumor  = ch_abra2_rna_input.tumor.ifEmpty(Channel.of([ [], [], [] ]))

        ABRA2_RNA(
            abra2_rna_normal,
            abra2_rna_tumor,
            fasta,
            gtf,
            junctions,
            []
        )
        // Collect versions, mix with channel outside scope.
        ch_versions = ch_versions.mix(ABRA2_RNA.out.versions)
        abra2_rna_out = ABRA2_RNA.out.abra2_normal.mix(ABRA2_RNA.out.abra2_tumor)
        ch_abra2_rna = ch_abra2_rna.mix(abra2_rna_out)

        // Allow user to save as CRAM
        if(!params.save_output_as_bam) {
            BAM_TO_CRAM_RNA(
                ch_abra2_rna,
                fasta,
                fasta_fai
            )
            ch_versions = ch_versions.mix(BAM_TO_CRAM_RNA.out.versions)
        }
    }

    emit:
        bam      = ch_abra2_rna
        cram     = ch_abra2_dna
        versions = ch_versions // channel: [ versions.yml ]
}
