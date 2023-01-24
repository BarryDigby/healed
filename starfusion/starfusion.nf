#!/usr/bin/env Nextflow

process starfusion {
// require:
//   FQS
//   params.starfusion$ctat_ref
//   params.starfusion$starfusion_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'starfusion_container'
  label 'starfusion'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/starfusion"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  path ctat_ref
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}*fusion_predictions.abridged.coding_effect*"), emit: coding_effect_fusions
  tuple val(pat_name), val(run), val(dataset), path("${dataset}*fusion_predictions.tsv"), emit: full_fusions

  script:
  """
  /usr/local/src/STAR-Fusion/STAR-Fusion --CPU ${task.cpus} --genome_lib_dir /\${PWD}/${ctat_ref} ${parstr} --left_fq ${fq1} --right_fq ${fq2} --output_dir .

  for i in `ls *fusion_predictions.abridged.coding_effect*`; do
    mv \${i} ${dataset}-${pat_name}-${run}.\${i}
  done

  for i in `ls *fusion_predictions.tsv`; do
    mv \${i} ${dataset}-${pat_name}-${run}.\${i}
  done
  """
}


process jstarfusion {
// require:
//   JUNCTIONS
//   params.starfusion$ctat_ref
//   params.starfusion$starfusion_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'starfusion_container'
  label 'starfusion'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/starfusion"

  input:
  tuple val(pat_name), val(run), val(dataset), path(junctions)
  path ctat_ref
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}*fusion_predictions.abridged.coding_effect*"), emit: coding_effect_fusions
  tuple val(pat_name), val(run), val(dataset), path("${dataset}*fusion_predictions.tsv"), emit: full_fusions

  script:
  """
  /usr/local/src/STAR-Fusion/STAR-Fusion --CPU ${task.cpus} --genome_lib_dir \${PWD}/${ctat_ref} ${parstr} -J ${junctions}  --output_dir .

  for i in `ls *fusion_predictions.abridged.coding_effect*`; do
    mv \${i} ${dataset}-${pat_name}-${run}.\${i}
  done

  for i in `ls *fusion_predictions.tsv`; do
    mv \${i} ${dataset}-${pat_name}-${run}.\${i}
  done
  """
}
