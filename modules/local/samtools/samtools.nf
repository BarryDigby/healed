#!/usr/bin/env nextflow

include { combine_sample_files } from '../utilities/utilities.nf'

process samtools_stats {
// Runs samtools stats
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(sam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.stats') - Output Alignment Stats File

// require:
//   ALNS
//   params.samtools$samtools_stats_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_stats'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/samtools_stats"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.stats'), emit: bam_stats

  script:
  """
  samtools stats ${parstr} ${aln} > ${aln}.stats
  """
}


process samtools_view {
// Runs samtools view
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(sam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.bam') - Output Alignment File

// require:
//   ALNS
//   params.samtools_view_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_view'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_view"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.bam'), emit: bams

  script:
  """
  aln_name=`echo ${aln}`
  bam_name=\${aln_name%.sam}.bam

  samtools view ${parstr} ${aln} > \${bam_name}
  """
}




process samtools_faidx {
// Runs samtools faidx
//
// input:
//   path fa -  Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: faidx_file
//     path(fa) -  Reference FASTA
//     path("${fa}.fai") - FASTA Index File

// require:
//   FA
//   params.samtools$samtools_faidx_parameters

  tag "${fa}"
  label 'samtools_container'
  label 'samtools_faidx'

  input:
  path(fa)
  val(parstr)

  output:
  tuple path("${fa}"), path("${fa}.fai"), emit: faidx_file

  script:
  """
  samtools faidx ${parstr} ${fa}
  """
}


process samtools_index {
// Runs samtools index
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(sam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.bam') - Output Alignment File

// require:
//   ALNS
//   params.samtools_index_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_index'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_index"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path(aln), path('*.bai'), emit: bams_and_bais
  tuple val(pat_name), val(run), val(dataset), path('*.bai'), emit: bais

  script:
  """
  samtools index ${parstr} ${aln} -@ ${task.cpus}
  """
}


process samtools_sort {
// Runs samtools sort
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(sam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.bam') - Output Alignment File

// require:
//   ALNS
//   params.samtools_sort_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_sort'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_sort"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.sorted.bam'), emit: bams

  script:
  """
  samtools sort ${aln} -o ${dataset}-${pat_name}-${run}.sorted.bam -@ ${task.cpus}
  """
}



process samtools_rmdup {
//runs samtools rmdup
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(aln) - Alignment File
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.rmdup.bam') - Removed Dups Alignment File

// require:
//   ALNS
//   params.samtools_rmdup_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_rmdup'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_rmdup"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.rmdup.bam'), emit: bams

  script:
  """
  ALN_NAME=`echo ${aln}`
  samtools rmdup ${aln} \${ALN_NAME%.bam}.rmdup.bam
  """
}


workflow sams_to_sorted_bams_w_indices {
// require:
//   ALNS
  take:
    alns
  main:
    samtools_view(
      alns,
      params.samtools$sams_to_sorted_bams_w_indices$samtools_view_parameters)
    samtools_sort(
      samtools_view.out.bams,
      params.samtools$sams_to_sorted_bams_w_indices$samtools_sort_parameters)
    samtools_index(
      samtools_sort.out.bams,
      params.samtools$sams_to_sorted_bams_w_indices$samtools_index_parameters)
  emit:
    bams = samtools_sort.out.bams
    bais = samtools_index.out.bais
    // Need to include bams_bais here
}


workflow bams_to_sorted_bams_w_indices {
// require:
//   ALNS
  take:
    bams
    samtools_index_parameters
  main:
    samtools_index(
      bams,
      samtools_index_parameters)
    combine_sample_files(
      bams,
      samtools_index.out.bais)
  emit:
    bams = bams
    bais = samtools_index.out.bais
    bams_bais = combine_sample_files.out.combined_set
}

process samtools_faidx_fetch {
// Runs samtools faidx to make a subset FASTA from BED file
//
// input:
//   path fa -  Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: faidx_file
//     path(fa) -  Reference FASTA
//     path("${fa}.fai") - FASTA Index File

// require:
//   FA
//   params.samtools$samtools_faidx_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_faidx'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_faidx"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bed)
  path fa
  val(suffix)
  val(parstr)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.${suffix}.fa"), emit: fetched_fastas

  script:
  """
  touch ${dataset}-${pat_name}-${run}.${suffix}.fa
  sed 's/\t/:/' ${bed} | sed 's/\t/-/' > ${dataset}-${pat_name}-${run}.samt.bed

  

  if [ -s ${dataset}-${pat_name}-${run}.samt.bed ]; then
    samtools faidx -r ${dataset}-${pat_name}-${run}.samt.bed ${parstr} ${fa} > ${dataset}-${pat_name}-${run}.${suffix}.fa
  fi
  """
}

process samtools_coverage {
// Runs samtools stats
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(sam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.stats') - Output Alignment Stats File

// require:
//   ALNS
//   params.samtools$samtools_coverage_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_coverage'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_coverage"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.coverage'), emit: bam_coverage

  script:
  """
  samtools coverage ${parstr} ${aln} > ${aln}.coverage
  """
}

process samtools_faidx_fetch_somatic {
// Runs samtools faidx to make a subset FASTA from BED file
//
// input:
//   path fa -  Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: faidx_file
//     path(fa) -  Reference FASTA
//     path("${fa}.fai") - FASTA Index File

// require:
//   FA
//   params.samtools$samtools_faidx_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'samtools_container'
  label 'samtools_faidx'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/samtools_faidx"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(beds)
  path fa
  val(suffix)
  val(parstr)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.${suffix}.fa"), emit: fetched_fastas

  script:
  """
  sed 's/\t/:/' ${bed} | sed 's/\t/-/' > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.samt.bed
  samtools faidx -r ${dataset}-${pat_name}-${norm_run}_${tumor_run}.samt.bed ${parstr} ${fa} > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.${suffix}.fa
  """
}

process samtools_faidx_fetch_somatic_folder {
// Runs samtools faidx to make a subset FASTA from BED file
//
// input:
//   path fa -  Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: faidx_file
//     path(fa) -  Reference FASTA
//     path("${fa}.fai") - FASTA Index File

// require:
//   FA
//   params.samtools$samtools_faidx_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'samtools_container'
  label 'samtools_faidx'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/samtools_faidx"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(beds)
  path fa
  val(suffix)
  val(parstr)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.exon_fas"), emit: fetched_fastas

  script:
  """
  mkdir -p ${dataset}-${pat_name}-${norm_run}_${tumor_run}.exon_fas
  mkdir -p tmp
  for i in `ls ${beds}`; do
    echo \${i}
    sed 's/\t/:/' ${beds}/\${i} | sed 's/\t/-/' > tmp/\${i%.bed}.samt.bed
    echo "Made tmp/\${i%.bed}.samt.bed"
    if [ -s tmp/\${i%.bed}.samt.bed ]; then
      samtools faidx -r tmp/\${i%.bed}.samt.bed ${parstr} ${fa} > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.exon_fas/\${i%.bed}.exons.fa
    fi
    echo "Made ${dataset}-${pat_name}-${norm_run}_${tumor_run}.exon_fas/\${i%.bed}.exons.fa"
  done
  echo "Done!"
  """
}


process samtools_faidx_fetch_folder {
// Runs samtools faidx to make a subset FASTA from BED file
//
// input:
//   path fa -  Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: faidx_file
//     path(fa) -  Reference FASTA
//     path("${fa}.fai") - FASTA Index File

// require:
//   FA
//   params.samtools$samtools_faidx_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'samtools_container'
  label 'samtools_faidx'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/samtools_faidx"

  input:
  tuple val(pat_name), val(run), val(dataset), path(beds)
  path fa
  val(suffix)
  val(parstr)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.exon_fas"), emit: fetched_fastas

  script:
  """
  mkdir -p ${dataset}-${pat_name}-${run}.exon_fas
  mkdir -p tmp
  for i in `ls ${beds}`; do
    echo \${i}
    sed 's/\t/:/' ${beds}/\${i} | sed 's/\t/-/' > tmp/\${i%.bed}.samt.bed
    echo "Made tmp/\${i%.bed}.samt.bed"
    if [ -s tmp/\${i%.bed}.samt.bed ]; then
      samtools faidx -r tmp/\${i%.bed}.samt.bed ${parstr} ${fa} > ${dataset}-${pat_name}-${run}.exon_fas/\${i%.bed}.exons.fa
    fi
    echo "Made ${dataset}-${pat_name}-${run}.exon_fas/\${i%.bed}.exons.fa"
  done
  echo "Done!"
  """
}


process bam_subsetter {

    tag "${dataset}/${pat_name}"

    label 'samtools_container'
    label 'samtools_view'

    publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/bam_subsetter"

    input:
    tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai), val(rna_run), path(rna_bam), path(rna_bai), path(bed)

    output:
    tuple val(pat_name), val(dataset), path("*nd*subset.bam"), path("*nd*subset.bam.bai"), path("*ad*subset.bam"), path("*ad*subset.bam.bai"), path("*ar*subset.bam"), path("*ar*subset.bam.bai"), emit: subsetted_bams
    tuple val(pat_name), val(dataset), val(norm_run), path("*nd*subset.bam"), path("*nd*subset.bam.bai"), val(tumor_run), path("*ad*subset.bam"), path("*ad*subset.bam.bai"), val(rna_run), path("*ar*subset.bam"), path("*ar*subset.bam.bai"), path(bed), emit: subsetted_bams_w_bed

    script:
    """
    for bam in `ls *bam | grep -v subset`; do
      samtools view \${bam} -bS -L ${bed} > \${bam%.bam}.subset.bam
      samtools index \${bam%.bam}.subset.bam
    done
    """
}
