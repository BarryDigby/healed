nextflow.enable.dsl=2

include { hlaprofiler_predict } from '../hlaprofiler.nf'

params.pwd = ''

workflow test_hlaprofiler_predict {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_hlaprofiler_predict/full.HLA00001N_HLA00002.A_1.fq.gz",
                "${params.pwd}/data/test_hlaprofiler_predict/full.HLA00001N_HLA00002.A_2.fq.gz"])
           .set{ test_hlaprofiler_predict_inputs }
    hlaprofiler_predict(
      test_hlaprofiler_predict_inputs,
      '')
}

workflow {
  test_hlaprofiler_predict()
}
