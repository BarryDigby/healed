include { FASTQC           } from '../../../modules/nf-core/fastqc/main'
include { FASTP            } from '../../../modules/nf-core/fastp/main'

workflow FASTQ_FASTQC_FASTP {
    take:
    reads            // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()

    if (!params.skip_fastqc) {
        FASTQC ( reads ).html.set { fastqc_html }
        ch_reports = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    if (params.trim_fastq || params.split_fastq > 0) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(reads,
                    [], // we are not using any adapter fastas at the moment
                    save_trimmed_fail,
                    save_merged)

            ch_reports = ch_reports.mix(
                                    FASTP.out.json.collect{meta, json -> json},
                                    FASTP.out.html.collect{meta, html -> html}
                                    )

            if (params.split_fastq){
                ch_reads_to_map = FASTP.out.reads.map{ key, reads ->

                        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                        [[
                            data_type:key.data_type,
                            id:key.id,
                            numLanes:key.numLanes,
                            patient: key.patient,
                            read_group:key.read_group,
                            sample:key.sample,
                            size:read_files.size(),
                            status:key.status,
                        ],
                        read_files]
                    }.transpose()
            }else{
                ch_reads_to_map = FASTP.out.reads
            }

            ch_versions = ch_versions.mix(FASTP.out.versions)
        } else {
            ch_reads_to_map = reads
        }

    emit:
    reads = ch_reads_to_map // channel: [ val(meta), [ reads ] ]
    reports = ch_reports.ifEmpty(null) // channel [zips, json, htmls.. ]
    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
