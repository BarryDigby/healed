//
// RNA MAPPING
//

include { CAT_FASTQ                                    } from '../../../modules/nf-core/cat/fastq/main'
include { STAR_ALIGN as STAR_ALIGN_SALMON              } from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS as SORT_STATS_SALMON } from '../../nf-core/bam_sort_stats_samtools'
include { SAMTOOLS_CONVERT as CONVERT_SALMON_BAM       } from '../../../modules/nf-core/samtools/convert/main'
include { STAR_ALIGN as STAR_ALIGN_STARFUSION          } from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS as SORT_STATS_STARFUSION } from '../../nf-core/bam_sort_stats_samtools'
include { SAMTOOLS_CONVERT as CONVERT_STARFUSION_BAM   } from '../../../modules/nf-core/samtools/convert/main'
include { STAR_ALIGN as STAR_ALIGN_ARRIBA              } from '../../../modules/nf-core/star/align/main'
include { BAM_SORT_STATS_SAMTOOLS as SORT_STATS_ARRIBA } from '../../nf-core/bam_sort_stats_samtools'
include { SAMTOOLS_CONVERT as CONVERT_ARRIBA_BAM       } from '../../../modules/nf-core/samtools/convert/main'


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
    // STAR FOR SALMON
    //
    ch_salmon_sorted_bam_index = Channel.empty()
    if('star_salmon' in rna_tools) {
        STAR_ALIGN_SALMON (
            reads,
            ch_map_index,
            gtf,
            false,
            params.seq_platform ?: '',
            'ILLUMINA'
        )
        ch_salmon_bam_transcript = STAR_ALIGN_SALMON.out.bam_transcript
        ch_versions              = ch_versions.mix(STAR_ALIGN_SALMON.out.versions.first())
        // Aligned.toTranscriptome.out.bam is never sorted.
        SORT_STATS_SALMON(ch_salmon_bam_transcript, fasta)
        ch_versions = ch_versions.mix(SORT_STATS_SALMON.out.versions)
        ch_reports  = ch_reports.mix(SORT_STATS_SALMON.out.stats)
        ch_reports  = ch_reports.mix(SORT_STATS_SALMON.out.flagstat)
        ch_reports  = ch_reports.mix(SORT_STATS_SALMON.out.idxstats)
        ch_salmon_sorted_bam_index = SORT_STATS_SALMON.out.bam.join(SORT_STATS_SALMON.out.bai)
        // Allow user to save as CRAM if they want.
        if(!params.save_output_as_bam && params.save_mapped) {
            CONVERT_SALMON_BAM(
                ch_salmon_sorted_bam_index,
                fasta,
                fasta_fai
            )
            ch_versions = ch_versions.mix(CONVERT_SALMON_BAM.out.versions)
        }
    }

    //
    // STAR FOR STAR-FUSION
    //
    ch_starfusion_sorted_bam_index = Channel.empty()
    ch_starfusion_junctions = Channel.empty()
    if('star_fusion' in rna_tools) {
        STAR_ALIGN_STARFUSION(
            reads,
            ch_map_index,
            gtf,
            false,
            params.seq_platform ?: '',
            'ILLUMINA'
        )
        ch_starfusion_bam       = STAR_ALIGN_STARFUSION.out.bam
        ch_starfusion_junctions = STAR_ALIGN_STARFUSION.out.junction
        ch_versions             = ch_versions.mix(STAR_ALIGN_STARFUSION.out.versions)
        // Sort downstream to optimise STAR alignment
        SORT_STATS_STARFUSION(ch_starfusion_bam, fasta)
        ch_versions = ch_versions.mix(SORT_STATS_STARFUSION.out.versions)
        ch_reports  = ch_reports.mix(SORT_STATS_STARFUSION.out.stats)
        ch_reports  = ch_reports.mix(SORT_STATS_STARFUSION.out.flagstat)
        ch_reports  = ch_reports.mix(SORT_STATS_STARFUSION.out.idxstats)
        ch_starfusion_sorted_bam_index = SORT_STATS_STARFUSION.out.bam.join(SORT_STATS_STARFUSION.out.bai)
        // Allow user to save as CRAM if they want.
        if(!params.save_output_as_bam && params.save_mapped) {
            CONVERT_STARFUSION_BAM(
                ch_starfusion_sorted_bam_index,
                fasta,
                fasta_fai
            )
            ch_versions = ch_versions.mix(CONVERT_STARFUSION_BAM.out.versions)
        }
    }

    //
    // STAR FOR ARRIBA
    //
    ch_arriba_sorted_bam_index = Channel.empty()
    if('arriba' in rna_tools) {
        STAR_ALIGN_ARRIBA(
            reads,
            ch_map_index,
            gtf,
            false,
            params.seq_platform ?: '',
            'ILLUMINA'
        )
        ch_arriba_bam = STAR_ALIGN_ARRIBA.out.bam
        ch_versions   = ch_versions.mix(STAR_ALIGN_ARRIBA.out.versions)
        // Sort downstream to optimise STAR alignment
        SORT_STATS_ARRIBA(ch_arriba_bam, fasta)
        ch_versions = ch_versions.mix(SORT_STATS_ARRIBA.out.versions)
        ch_reports  = ch_reports.mix(SORT_STATS_ARRIBA.out.stats)
        ch_reports  = ch_reports.mix(SORT_STATS_ARRIBA.out.flagstat)
        ch_reports  = ch_reports.mix(SORT_STATS_ARRIBA.out.idxstats)
        ch_arriba_sorted_bam_index = SORT_STATS_ARRIBA.out.bam.join(SORT_STATS_ARRIBA.out.bai)
        // Allow user to save as CRAM if they want.
        if(!params.save_output_as_bam && params.save_mapped) {
            CONVERT_ARRIBA_BAM(
                ch_arriba_sorted_bam_index,
                fasta,
                fasta_fai
            )
            ch_versions = ch_versions.mix(CONVERT_ARRIBA_BAM.out.versions)
        }
    }

    emit:
        arriba_bam            = ch_arriba_sorted_bam_index // tuple meta, bam, bai [sorted (not required but should improve compute) bam]
        salmon_bam_transcript = ch_salmon_sorted_bam_index // tuple meta, bam, bai [sorted transcriptome bam file]
        starfusion_bam        = ch_starfusion_sorted_bam_index // tuple meta(val), bam, bai [sorted bam]
        starfusion_junctions  = ch_starfusion_junctions    // tuple meta, out.tab
        versions              = ch_versions
        reports               = ch_reports
}
