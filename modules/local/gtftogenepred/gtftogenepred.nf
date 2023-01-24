#!/usr/bin/env nextflow

process gtf_to_genepred {
// require:
//   GTF
//   params.gtf_to_genepred$gtf_to_genepred_parameters

  tag "${gtf}"
  label 'gtftogenepred_container'
  label 'gtftogenepred'

  input:
  path gtf
  val parstr

  output:
  path "*genepred", emit: genepred

  script:
  """
  AGTF=\$(echo ${gtf})
  gtfToGenePred ${parstr} ${gtf} \${AGTF%.gtf}.genepred
  """
}
