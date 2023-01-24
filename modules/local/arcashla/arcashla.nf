#!/usr/bin/env nextflow

process arcashla_extract {

  tag "${dataset}/${pat_name}/${prefix}"
  label 'arcashla_container'
  label 'arcashla_extract'

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(aln), path(idx)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("*1.fq.gz"), path("*2.fq.gz"), emit: extd_fqs

  script:
  """
  arcasHLA extract ${parstr} ${aln}
  """
}


process arcashla_genotype {

  tag "${dataset}/${pat_name}/${prefix}"
  label 'arcashla_container'
  label 'arcashla_genotype'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${prefix}/arcashla_genotype"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("*.genotype.json), emit: alleles

  script:
  """
  arcasHLA genotype ${parstr} ${fq1} ${fq2}
  """
}
