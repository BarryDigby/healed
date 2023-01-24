nextflow.enable.dsl=2

include { abra2 } from '../abra2.nf'
include { abra2_cadabra } from '../abra2.nf'

params.pwd = ''

workflow test_abra2 {
  main:
    Channel.of(['test_patient',
                'test_dataset',
                'test_tumor_run',
                "${params.pwd}/data/test_abra2/test.n.bam",
                "${params.pwd}/data/test_abra2/test.n.bam.bai",
                'test_norm_run',
                "${params.pwd}/data/test_abra2/test.t.bam",
                "${params.pwd}/data/test_abra2/test.t.bam.bai"])
           .set{ sample_inputs }
    Channel.of(["${params.pwd}/data/test_abra2/chr1.fa",
                "${params.pwd}/data/test_abra2/dummy_file"])
           .set{ ref_inputs }
    abra2(
      sample_inputs,
      ref_inputs,
      "${params.pwd}/data/test_abra2/test.bed",
      '',
      '.abra2')
}

workflow test_abra2_cadabra {
  main:
    Channel.of(['test_patient',
                'test_dataset',
                'test_tumor_run',
                "${params.pwd}/data/test_abra2_cadabra/test.n.norm_abra.bam",
                "${params.pwd}/data/test_abra2_cadabra/test.n.norm_abra.bam.bai",
                'test_norm_run',
                "${params.pwd}/data/test_abra2_cadabra/test.t.tumor_abra.bam",
                "${params.pwd}/data/test_abra2_cadabra/test.t.tumor_abra.bam.bai"])
           .set{ sample_inputs }
    abra2_cadabra(
      sample_inputs,
      "${params.pwd}/data/test_abra2/chr1.fa",
      '.abra2',
      '')
}
