#!/usr/bin/env nextflow

process star_index {
// Runs star --runMode genomeGenerate
//
// input:
//   path fa - Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: idx_files
//     path(fa) - Reference FASTA
//     path("index_files/*") - Index Files

// require:
//   params.star$star_reference
//   params.star$star_index_parameters

  tag "${fa}"
  label 'star_container'
  label 'star_index'
  storeDir "${params.shared_dir}/${fa}/star_index"

  input:
  path fa
  val parstr

  output:
  tuple path(fa), path("index_files/*"), emit: idx_files

  script:
  """
  echo "${params.shared_dir}/${fa}/star_index"
  mkdir -p index_files
  STAR \
    --runMode genomeGenerate \
    --genomeDir index_files \
    --genomeFastaFiles ${fa} \
    ${parstr} \
    --runThreadN ${task.cpus}
  """
}


process star_map {
// Runs STAR in mapping mode
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: alns
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*am') - Alignment File

// require:
//   FQS
//   IDX_FILES
//   params.star$star_map_parameters
//   params.star$star_alt_capture

  tag "${dataset}/${pat_name}/${run}"
  label 'star_container'
  label 'star_map'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/star_map"
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/star_map", pattern: "*Log.final.out"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx)
  val parstr
  val alt_capture
  path gtf

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*am"), emit: alns
  tuple val(pat_name), val(run), val(dataset), path("SJ.out.tab"), emit: standard_junctions
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}.Log.final.out"), emit: star_logs
  tuple val(pat_name), val(run), val(dataset), path("alt_cap_alns/${dataset}-${pat_name}-${run}*am"), emit: alt_alns optional true
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}.Chimeric.out.junction"), emit: junctions optional true

  script:
  def gtf_proxy  = gtf.name != 'dummy_file' ? "--sjdbGTFfile ${gtf}" : ''
  """
  STAR \
    --genomeDir . \
    --readFilesCommand zcat \
    --readFilesIn ${fq1} ${fq2} \
    ${gtf_proxy} \
    ${parstr} \
    --limitBAMsortRAM ${task.memory.toBytes()} \
    --runThreadN ${task.cpus} \
    --outSAMattrRGline ID:${dataset}-${pat_name}-${run} SM:${dataset}-${pat_name}-${run}

  for i in \$(ls *junction); do
      mv \${i} ${dataset}-${pat_name}-${run}.Chimeric.out.junction
  done

  for i in \$(ls *am); do
      ALIGN=\${i}
      SUFFIX=\$(echo \${ALIGN} | cut -f 1- -d .)
      mv \${ALIGN} ${dataset}-${pat_name}-${run}.\${SUFFIX}
  done

  ALT_CAPTURE_ECHO="${alt_capture}"

  if [[ ! -z \${ALT_CAPTURE_ECHO} ]]; then
    #This is a dumb hack, but should work.
    mkdir -p alt_cap_alns
    for aln in `ls *am | grep "${alt_capture}"`; do
      mv \${aln} alt_cap_alns/
    done
  fi

  mv Log.final.out ${dataset}-${pat_name}-${run}.Log.final.out
  """
}
