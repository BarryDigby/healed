#!/usr/bin/env nextflow

process multiqc {
// require:
//   METRIC_FILES
//   params.multiqc$multiqc_parameterss
  
  label 'multiqc_container'
  label 'multiqc'
  publishDir "${params.qc_out_dir}/multiqc"

  input:
  path metric_files
  val parstr

  output:
  path "multiqc*", emit: multiqc_reports
  path "multiqc_data", emit: multiqc_data

  script:
  """
  multiqc . ${parstr}
  """
}
