nextflow.enable.dsl=2

include { netmhcpan } from '../netmhcpan.nf'
include { netmhcpan_rna } from '../netmhcpan.nf'

params.pwd = ''

workflow test_netmhcpan {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tumor_run',
                'test_dataset',
                "${params.pwd}/data/test_netmhcpan/test.fa",
                "${params.pwd}/data/test_netmhcpan/test.alleles"])
           .set{ sample_inputs }
    netmhcpan(
      sample_inputs,
      '')
}

