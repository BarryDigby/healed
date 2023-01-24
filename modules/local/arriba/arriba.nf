#!/usr/bin/env Nextflow

process arriba {
// require:
//   BAM
//   params.arriba$arriba_ref
//   params.starfusion$arriba_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'arriba_container'
  label 'arriba'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/arriba"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)
  path arriba_ref
  path fa
  path gtf
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*fusions.tsv"), emit: fusions

  script:
  """
  /arriba*/arriba \
    -x ${bam} \
    -o ${dataset}-${pat_name}-${run}.arriba_fusions.tsv \
    -O ${dataset}-${pat_name}-${run}.arriba_fusions.discarded.tsv \
    -a ${fa} \
    -g ${gtf} \
    -b ${arriba_ref}/blacklist_hg38_GRCh38_v2.3.0.tsv.gz \
    -k ${arriba_ref}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz \
    -t ${arriba_ref}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz \
    -p ${arriba_ref}/protein_domains_hg38_GRCh38_v2.3.0.gff3
  """
}

process arriba_star_map {
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
  label 'arriba_container'
  label 'arriba_star_map'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx)
  val parstr
  path gtf

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*am"), emit: alns

  script:
  """
  STAR \
    --genomeDir . \
    --readFilesCommand zcat \
    --readFilesIn ${fq1} ${fq2} \
    --sjdbGTFfile ${gtf} \
    ${parstr} \
    --limitBAMsortRAM ${task.memory.toBytes()} \
    --runThreadN ${task.cpus} \
    --outSAMattrRGline ID:${dataset}-${pat_name}-${run} SM:${dataset}-${pat_name}-${run}

  mv Aligned.sortedByCoord.out.bam ${dataset}-${pat_name}-${run}.Aligned.sortedByCoord.out.bam
  """
}

process arriba_star_index {
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
  label 'arriba_container'
  label 'arriba_star_index'
  storeDir "${params.shared_dir}/${fa}/star_index/arriba"

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
