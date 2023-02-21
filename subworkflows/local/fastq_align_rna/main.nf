//
// RNA MAPPING
//

include { CAT_FASTQ                      } from '../../../modules/nf-core/cat/fastq/main'
include { STAR_ALIGN as STAR_ALIGN_QUANT } from '../../../modules/nf-core/star/align/main'

workflow FASTQ_ALIGN_RNA {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: path star/
        gtf          // channel: [mandatory] gtf
        rna_tools    // array list

    main:

    ch_versions = Channel.empty()

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
    // Stage empty channels, reduce these as you see fit
    //
    ch_orig_bam       = Channel.empty()
    ch_log_final      = Channel.empty()
    ch_log_out        = Channel.empty()
    ch_log_progress   = Channel.empty()
    ch_bam_sorted     = Channel.empty()
    ch_bam_transcript = Channel.empty()
    ch_fastq          = Channel.empty()
    ch_tab            = Channel.empty()
    if('star_salmon' in rna_tools) {
        STAR_ALIGN_QUANT (
            reads,
            ch_map_index,
            gtf,
            false,
            params.seq_platform ?: '',
            'ILLUMINA'
        )
        ch_orig_bam       = STAR_ALIGN_QUANT.out.bam
        ch_log_final      = STAR_ALIGN_QUANT.out.log_final
        ch_log_out        = STAR_ALIGN_QUANT.out.log_out
        ch_log_progress   = STAR_ALIGN_QUANT.out.log_progress
        ch_bam_sorted     = STAR_ALIGN_QUANT.out.bam_sorted
        ch_bam_transcript = STAR_ALIGN_QUANT.out.bam_transcript
        ch_fastq          = STAR_ALIGN_QUANT.out.fastq
        ch_tab            = STAR_ALIGN_QUANT.out.tab
        ch_versions       = ch_versions.mix(STAR_ALIGN_QUANT.out.versions.first())
    }

    emit:
        bam       = ch_orig_bam // [channel] tuple
        junctions = ch_tab  // [channel] tuple
        versions  = ch_versions   // [channel] path
}
