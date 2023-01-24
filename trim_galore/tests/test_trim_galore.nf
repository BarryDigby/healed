nextflow.enable.dsl=2

include { trim_galore } from '../trim_galore.nf'

params.pwd = ''

workflow test_trim_galore {
  main:
    Channel.of(['test_patient', 'test_run', 'test_dataset', "${params.pwd}/data/test_37.5k_1.fastq.gz", "${params.pwd}/data/test_37.5k_2.fastq.gz"]).set{ trim_galore_inputs }
    trim_galore(
      trim_galore_inputs,
      '--paired -a A{10}')
}

workflow {
  test_trim_galore()
}
