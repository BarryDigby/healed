nextflow.enable.dsl=2

include { netmhcstabpan } from '../netmhcstabpan.nf'
include { netmhcstabpan_rna } from '../netmhcstabpan.nf'

params.pwd = ''

workflow test_netmhcstabpan {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tumor_run',
                'test_dataset',
                "${params.pwd}/data/test_netmhcstabpan/test.fa",
                "${params.pwd}/data/test_netmhcstabpan/test.alleles"])
           .set{ sample_inputs }
    netmhcstabpan(
      sample_inputs,
      '')
}

