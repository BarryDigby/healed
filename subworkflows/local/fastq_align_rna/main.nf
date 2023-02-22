//
// RNA MAPPING
//

include { CAT_FASTQ                                   }    from '../../../modules/nf-core/cat/fastq/main'
include { STAR_ALIGN as STAR_ALIGN_QUANT              }    from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS as SORT_STATS_QUANT }    from '../../nf-core/bam_sort_stats_samtools'
include { SAMTOOLS_CONVERT as CONVERT_QUANT           }    from '../../../modules/nf-core/samtools/convert/main'
include { STAR_ALIGN as STAR_ALIGN_FUSION             }    from '../../../modules/nf-core/star/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as INDEX_STATS_FUSION } from '../bam_index_stats_samtools'
include { SAMTOOLS_CONVERT as CONVERT_FUSION          }    from '../../../modules/nf-core/samtools/convert/main'

workflow FASTQ_ALIGN_RNA {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: path star/
        fasta        // channel: path fasta
        fasta_fai    // channel: path fai
        gtf          // channel: [mandatory] gtf
        rna_tools    // array list

    main:

    ch_versions = Channel.empty()
    ch_reports = Channel.empty()

    // Do not split the FASTQ files for RNA seq, the sheer volume of output
    // files from STAR makes this a headache. Recall that RNASeq did away
    // with the concept of lanes, instead opting to provide duplicate sample
    // id (first column) to indicate the presence of 'lanes'.
    //
    // BUT, we are adopting a lot of sarek logic. Consider the below channel
    // whereby a sample is indeed split over two lanes as indicated by our
    // samplesheet:
    //
    // [[patient:GM12878, assay:rna, status:normal, sample:GM12878N_T1, lane:1, id:GM12878N-1, numLanes:1, read_group:"@RG\tID:null.GM12878N_T1.1\tPU:1\tSM:GM12878_GM12878N_T1\tLB:GM12878N_T1\tDS:dev_data/chr21.fa\tPL:ILLUMINA", data_type:fastq, size:1, strandedness:reverse], [/data/github/healed/dev_data/GM12878_R1.fastq.gz, /data/github/healed/dev_data/GM12878_R2.fastq.gz]]
    // [[patient:GM12878, assay:rna, status:normal, sample:GM12878N_T2, lane:2, id:GM12878N-2, numLanes:1, read_group:"@RG\tID:null.GM12878N_T2.2\tPU:2\tSM:GM12878_GM12878N_T2\tLB:GM12878N_T2\tDS:dev_data/chr21.fa\tPL:ILLUMINA", data_type:fastq, size:1, strandedness:reverse], [/data/github/healed/dev_data/GM128781_R1.fastq.gz, /data/github/healed/dev_data/GM128781_R2.fastq.gz]]
    //
    // compared to just one sample on one lane:
    //
    // [[patient:GM12878, assay:rna, status:normal, sample:GM12878N_T1, lane:1, id:GM12878N-1, numLanes:1, read_group:"@RG\tID:null.GM12878N_T1.1\tPU:1\tSM:GM12878_GM12878N_T1\tLB:GM12878N_T1\tDS:dev_data/chr21.fa\tPL:ILLUMINA", data_type:fastq, size:1, strandedness:reverse], [/data/github/healed/dev_data/GM12878_R1.fastq.gz, /data/github/healed/dev_data/GM12878_R2.fastq.gz]]
    //
    // Gameplan here is to strip -1/-2 etc or _T1, _T2 to use groupTuple():

    ch_reads.map{ meta, reads ->

        new_id         = meta.id - ~/_(T\d+)-(\d+)/
        new_sample     = meta.sample - ~/_T\d+/

        // drop lane & size info + RG (STAR adds this for us)
        [[
            id:           new_id,
            sample:       new_sample,
            patient:      meta.patient,
            status:       meta.status,
            strandedness: meta.strandedness
            ],
        reads]
    }
    .groupTuple()
    .branch {
        meta, reads ->
            single  : reads.size() == 1
                return [ meta, reads.flatten() ]
            multiple: reads.size() > 1
                return [ meta, reads.flatten() ]
    }
    .set { ch_fastq }

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { reads }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // STAR FOR QUANTIFICATION
    //
    ch_quant_sorted_bam_index = Channel.empty()
    if('star_salmon' in rna_tools) {
        STAR_ALIGN_QUANT (
            reads,
            ch_map_index,
            gtf,
            false,
            params.seq_platform ?: '',
            'ILLUMINA'
        )
        ch_quant_bam_transcript = STAR_ALIGN_QUANT.out.bam_transcript
        ch_versions             = ch_versions.mix(STAR_ALIGN_QUANT.out.versions.first())
        // Aligned.toTranscriptome.out.bam is never sorted.
        SORT_STATS_QUANT(ch_quant_bam_transcript, fasta)
        ch_versions = ch_versions.mix(SORT_STATS_QUANT.out.versions)
        ch_reports  = ch_reports.mix(SORT_STATS_QUANT.out.stats)
        ch_reports  = ch_reports.mix(SORT_STATS_QUANT.out.flagstat)
        ch_reports  = ch_reports.mix(SORT_STATS_QUANT.out.idxstats)
        ch_quant_sorted_bam_index = SORT_STATS_QUANT.out.bam.join(SORT_STATS_QUANT.out.bai)
        // Allow user to save as CRAM if they want.
        if(!params.save_output_as_bam && params.save_mapped) {
            CONVERT_QUANT(
                ch_quant_sorted_bam_index,
                fasta,
                fasta_fai
            )
            ch_versions = ch_versions.mix(CONVERT_QUANT.out.versions)
        }
    }

    //
    // STAR FOR STAR-FUSION
    //
    ch_fusion_sorted_bam = Channel.empty()
    if('star_fusion' in rna_tools) {
        STAR_ALIGN_FUSION(
            reads,
            ch_map_index,
            gtf,
            false,
            params.seq_platform ?: '',
            'ILLUMINA'
        )
        ch_fusion_sorted_bam = STAR_ALIGN_FUSION.out.bam_sorted
        ch_versions          = ch_versions.mix(STAR_ALIGN_FUSION.out.versions)
        INDEX_STATS_FUSION(ch_fusion_sorted_bam, fasta)
        ch_versions = ch_versions.mix(INDEX_STATS_FUSION.out.versions)
        ch_reports  = ch_reports.mix(INDEX_STATS_FUSION.out.stats)
        ch_reports  = ch_reports.mix(INDEX_STATS_FUSION.out.flagstat)
        ch_reports  = ch_reports.mix(INDEX_STATS_FUSION.out.idxstats)
        ch_fusion_sorted_bam_index = INDEX_STATS_FUSION.out.bam.join(INDEX_STATS_FUSION.out.bai)
        // Allow user to save as CRAM if they want.
        if(!params.save_output_as_bam && params.save_mapped) {
            CONVERT_FUSION(
                ch_fusion_sorted_bam_index,
                fasta,
                fasta_fai
            )
            ch_versions = ch_versions.mix(CONVERT_FUSION.out.versions)
        }
    }


    emit:
        quant_bam_transcript = ch_quant_sorted_bam_index // tuple meta, bam, bai [sorted transcriptome bam file]
        versions             = ch_versions
        reports              = ch_reports
}
