nextflow.enable.dsl=2

include { picard_mark_duplicates } from '../picard2.nf'
include { picard_create_seq_dict } from '../picard2.nf'

params.pwd = ''

workflow test_picard_mark_duplicates {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_picard_mark_duplicates/test.t.bam"])
           .set{ sample_inputs }
    picard_mark_duplicates(
      sample_inputs,
      '')
}

workflow test_picard_create_seq_dict {
  main:
    picard_create_seq_dict(
      "${params.pwd}/data/test_picard_create_seq_dict/chr1.fa",
      '')
}

