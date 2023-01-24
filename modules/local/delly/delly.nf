#!/usr/bin/env nextflow

process delly_lr {
// Runs Delly for long reads
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
//   params.cutesv$cutesv_mem_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'delly_container'
  label 'delly_lr'
  cache 'lenient'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  path fa
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.vcf'), emit: vcfs

  script:
  """
  delly lr \
    ${parstr} \
    -g ${fa} \
    ${bam} >\
    ${dataset}-${pat_name}-${run}.delly.vcf
  """
}
