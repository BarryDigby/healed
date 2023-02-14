nextflow.enable.dsl=2

include { snpeff_ann } from '../snpeff.nf'
include { snpeff_ann_germline } from '../snpeff.nf'

params.pwd = ''

workflow test_snpeff_ann {
  main:
    Channel.of(['test_patient',
                'norm_test_run',
                'tumor_test_run',
                'test_dataset',
                "${params.pwd}/data/test_snpeff_ann/test.vcf.gz"])
           .set{ sample_inputs }
    Channel.of("${params.pwd}/data/test_snpeff_ann/GRCh38.GENCODEv37")
           .set{ ref_inputs }
    snpeff_ann(
      sample_inputs,
      ref_inputs)
}

workflow test_snpeff_ann_germline {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_snpeff_ann_germline/test.vcf.gz"])
           .set{ sample_inputs }
    Channel.of("${params.pwd}/data/test_snpeff_ann_germline/GRCh38.GENCODEv37")
           .set{ ref_inputs }
    snpeff_ann_germline(
      sample_inputs,
      ref_inputs)
}
