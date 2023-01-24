#!/usr/bin/env nextflow

process mixcr_align {
// Runs mixcr align
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: vdjca_files
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*vdjca') - VDJCA Output Files
//
// require:
//   FQS
//   params.mixcr$mixcr_align_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'mixcr_container'
  label 'mixcr_align'
  tuple val(pat_name), val(run), val(dataset), path('*vdjca'), emit: vdjca_files

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*vdjca'), emit: vdjca_files

  script:
  """
  mixcr align ${parstr} ${fq1} ${fq2} ${dataset}-${pat_name}-${run}.vdjca
  """
}


process mixcr_assemble {
// Runs mixcr assemble
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(vdjca) - VDJCA Input File
//   val if_config - null = assemble, non-null = assembleContigs
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: cln_files
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.cln*') - Output File (either *.clna or *.clns)
//
// require:
//   VDJCAS
//   params.mixcr$mixcr_assemble_if_config
//   params.mixcr$mixcr_assemble_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'mixcr_container'
  label 'mixcr_assemble'

  input:
  tuple val(pat_name), val(run), val(dataset), path(vdjca)
  val if_config
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.cln*'), emit: cln_files

  script:
  """
  SUB=''
  if [ -z ${if_config} ]; then
    SUB='assemble'
  else
    SUB='assembleContigs'
  fi

  SUFFIX='.clns'
  PARAMS=`echo ${parstr}`
  if [ ! -z \${PARAMS} ]; then
    if [[ \${PARAMS} =~ ([[:space:]]-a[[:space:]]) || \${PARAMS} =~ (--write-alignments) ]]; then
      SUFFIX='.clna'
    fi
  fi

  mixcr \${SUB} ${parstr} ${vdjca} ${dataset}-${pat_name}-${run}\${SUFFIX}
  """
}


process mixcr_export {
// Runs mixcr export
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(infile) - Input File (either *.clna or *.clns)
//   val exports - "clones" = exportClones, *.clna and "alignments" = exportAlignments
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: exported_files
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*txt') - Output file
//
// require:
//   INPUT_FILES
//   params.mixcr$mixcr_export_exports
//   params.mixcr$mixcr_export_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'mixcr_container'
  label 'mixcr_export'

  input:
  tuple val(pat_name), val(run), val(dataset), path(infile)
  val exports
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*txt'), emit: exported_files

  script:
  """
  LEXPORTS=`echo ${exports}`
  if [[ \${LEXPORTS} =~ (clones) || -z \${LEXPORTS} ]];then
    mixcr exportClones ${parstr} ${infile} ${dataset}-${pat_name}-${run}.clones.txt
  elif [[ ${infile} =~ (.clna) && (\${LEXPORTS} =~ alignments || -z \${LEXPORTS}) ]];then
    mixcr exportAlignments ${parstr} ${infile} ${dataset}-${pat_name}-${run}.alignments.txt
  fi
  """
}


process mixcr_shotgun {
// Runs mixcr analyze shotgun & exportReads
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(datasest) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ2
//   val regex - Output Regular Expression
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: clonotypes
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(datasest) - Dataset
//     path("${regex}") - exportClones Clonotypes
//   tuple => emit: report
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(datasest) - Dataset
//     path("*.report") - Combined Report
//   tuple => emit: aligned_r1
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(datasest) - Dataset
//     path("*_aligned_r1.fq.gz" - mixcr align R1 Reads
//   tuple => emit: aligned_r2
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(datasest) - Dataset
//     path("*_aligned_r2.fq.gz" - mixcr align R2 Reads
//
// require:
//   FQS
//   params.mixcr$mixcr_shotgun_regex
//   params.mixcr$mixcr_shotgun_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'mixcr_container'
  label 'mixcr_shotgun'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val regex
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${regex}"), emit: clonotypes             // exportClones clonotypes text file
  tuple val(pat_name), val(run), val(dataset), path("*.report"), emit: report                 // combined report from align, assemblePartial, extend, assemble and export
  tuple val(pat_name), val(run), val(dataset), path("*_aligned_r1.fq.gz"), emit: aligned_r1   // aligned r1 reads from mixcr align
  tuple val(pat_name), val(run), val(dataset), path("*_aligned_r2.fq.gz"), emit: aligned_r2   // aligned r2 reads from mixcr align

  script:
  """
  mixcr analyze shotgun -t ${task.cpus} --align "-OsaveOriginalReads=true" ${parstr} ${fq1} ${fq2} ${dataset}-${pat_name}-${run}; \
  mixcr exportReads ${dataset}-${pat_name}-${run}.vdjca ${dataset}-${pat_name}-${run}_aligned_r1.fq.gz ${dataset}-${pat_name}-${run}_aligned_r2.fq.gz
  """
}
