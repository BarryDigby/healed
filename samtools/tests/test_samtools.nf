nextflow.enable.dsl=2

include { samtools_coverage } from '../samtools.nf'
include { samtools_faidx } from '../samtools.nf'
include { samtools_index } from '../samtools.nf'
include { samtools_rmdup } from '../samtools.nf'
include { samtools_sort } from '../samtools.nf'
include { samtools_stats } from '../samtools.nf'
include { samtools_view } from '../samtools.nf'

params.pwd = ''

workflow test_samtools_coverage {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_samtools_coverage/test.bam"])
           .set{ sample_inputs }
    samtools_coverage(
      sample_inputs,
      '')
}

workflow test_samtools_index {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_samtools_index/test.bam"])
           .set{ sample_inputs }
    samtools_index(
      sample_inputs,
      '')
}

workflow test_samtools_faidx {
  main:
    samtools_faidx(
      "${params.pwd}/data/test_samtools_faidx/chr1.fa",
      '')
}

workflow test_samtools_rmdup {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_samtools_index/test.bam"])
           .set{ sample_inputs }
    samtools_rmdup(
      sample_inputs,
      '')
}

workflow test_samtools_sort {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_samtools_index/test.bam"])
           .set{ sample_inputs }
    samtools_sort(
      sample_inputs,
      '')
}

workflow test_samtools_stats {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_samtools_index/test.bam"])
           .set{ sample_inputs }
    samtools_stats(
      sample_inputs,
      '')
}

workflow test_samtools_view {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_samtools_index/test.bam"])
           .set{ sample_inputs }
    samtools_view(
      sample_inputs,
      '')
}
