nextflow.enable.dsl=2

include { pyclonevi_fit } from '../pyclone-vi.nf'

params.pwd = ''

workflow test_pyclonevi_fit {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tumor_run',
                'test_dataset',
                "${params.pwd}/data/test_pyclonevi_fit/test.pcvi_input"])
           .set{ sample_inputs }
    pyclonevi_fit(
      sample_inputs,
      '')
}
