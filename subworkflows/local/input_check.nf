//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def create_fastq_channel(LinkedHashMap row) {

    def meta = [:]
    meta.patient      = row.patient
    meta.assay        = row.assay
    meta.status       = row.status
    meta.sample       = row.sample
    meta.lane         = row.lane
    if (meta.assay == 'rna') { meta.strandedness = row.strandedness }
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    
    fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    
    return fastq_meta
}
