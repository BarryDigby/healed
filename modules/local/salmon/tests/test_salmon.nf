nextflow.enable.dsl=2

include { salmon_index } from '../salmon.nf'
include { salmon_map_quant } from '../salmon.nf'
include { salmon_aln_quant } from '../salmon.nf'

params.pwd = ''

workflow test_salmon_index{
  main:
    salmon_index(
      "${params.pwd}/data/test_salmon_index/test.tx.fa",
      '')
}

workflow test_salmon_map_quant {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_salmon_map_quant/test_37.5k_1.fastq.gz",
                "${params.pwd}/data/test_salmon_map_quant/test_37.5k_2.fastq.gz"])
           .set{ sample_inputs }
    Channel.of(["${params.pwd}/data/test_salmon_map_quant/test.tx.fa",
                "${params.pwd}/null/test.tx.fa/salmon_index/test.tx.fa.index"])
           .set{ ref_inputs }
    salmon_map_quant(
      sample_inputs,
      ref_inputs,
      '')
}

workflow test_salmon_aln_quant {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_salmon_aln_quant/test.bam"])
           .set{ sample_inputs }
    salmon_aln_quant(
      sample_inputs,
      "${params.pwd}/data/test_salmon_aln_quant/test.tx.fa",
      '')
}
