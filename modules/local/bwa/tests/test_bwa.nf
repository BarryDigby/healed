nextflow.enable.dsl=2

include { bwa_index } from '../bwa.nf'
include { bwa_mem_samtools_sort } from '../bwa.nf'

params.pwd = ''

workflow test_bwa_index {
  main:
    bwa_index(
      "${params.pwd}/data/test_bwa_index/NC_045512.2.fa",
      '')
}

workflow test_bwa_mem_samtools_sort {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_bwa_mem_samtools_sort/V7980_S98_R1_001.fastq.gz",
                "${params.pwd}/data/test_bwa_mem_samtools_sort/V7980_S98_R2_001.fastq.gz"])
           .set{ sample_inputs }
    Channel.of(["${params.pwd}/data/test_bwa_mem_samtools_sort/NC_045512.2",
                ["${params.pwd}/data/test_bwa_mem_samtools_sort/idx/NC_045512.2.amb",
                 "${params.pwd}/data/test_bwa_mem_samtools_sort/idx/NC_045512.2.ann",
                 "${params.pwd}/data/test_bwa_mem_samtools_sort/idx/NC_045512.2.bwt",
                 "${params.pwd}/data/test_bwa_mem_samtools_sort/idx/NC_045512.2.pac",
                 "${params.pwd}/data/test_bwa_mem_samtools_sort/idx/NC_045512.2.sa"]])
           .set{ ref_inputs }
    bwa_mem_samtools_sort(
      sample_inputs,
      ref_inputs,
      '')
}
