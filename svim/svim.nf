#!/usr/bin/env nextflow

process svim_alignment {
// Runs svim in alignment mode
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq1) - FASTQ 2
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//   parstr - Additional Parameters
//
// output:
//   tuple => emit: sams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.sam') - Output SAM File

// require:
//   FQS
//   IDX_FILES
//   params.svim$svim_mem_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'svim_container'
  label 'svim'
  cache 'lenient'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  path fa
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.vcf'), emit: vcfs

  script:
  """
  svim alignment \
    ${dataset}-${pat_name}-${run}.svim_outputs \
    ${bam} \
    ${fa} \
    ${parstr}

 cp ${dataset}-${pat_name}-${run}.svim_outputs/variants.vcf ${dataset}-${pat_name}-${run}.svim.vcf
  """
}
