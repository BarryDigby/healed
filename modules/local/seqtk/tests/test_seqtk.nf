nextflow.enable.dsl=2

include { seqtk_sample } from '../seqtk.nf'

params.pwd = ''

workflow test_seqtk_sample {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_seqtk_sample/test_1.fastq.gz",
                "${params.pwd}/data/test_seqtk_sample/test_2.fastq.gz"])
           .set{ sample_inputs }
    seqtk_sample(
      sample_inputs,
      100000,
      '1234',
      '.subd',
      '')
}
