nextflow.enable.dsl=2

include { gffread_make_tx_fa } from '../gffread.nf'

params.pwd = ''

workflow test_gffread_make_tx_fa {
  main:
    gffread_make_tx_fa(
      "${params.pwd}/data/test_gffread_make_tx_fa/test.fa",
      "${params.pwd}/data/test_gffread_make_tx_fa/test.gtf")
}

workflow {
  test_gffread_make_tx_fa()
}
