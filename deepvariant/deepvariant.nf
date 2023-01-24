#!/usr/bin/env nextflow

process deepvariant {
// Runs DeepVariant on a sample
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     path(bam) - BAM File
//     path(bai) - BAI File
//   tuple
//     path(fa) - Reference FASTA
//     path(fai) - Reference FASTA Index
//     path(dict) - Sequence Dictionary File
//   path bed - BED File
//   val model_type - Model Type
//   val params - Additional Parameters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: germline_vcfs
//     val(pat_name) - Patient Name
//     val(run_name) - Run Name
//     val(dataset) - Dataset
//     path("*germline.vcf") - VCF

// require:
//   BAMS
//   REFS_W_INDICES
//   params.deepvariant$targets_bed
//   params.deepvariant$deepvariant_model_type
//   params.deepvariant$deepvariant_parameters
//   params.deepvariant$deepvariant_suffix

  tag "${dataset}/${pat_name}/${run}"
  label 'deepvariant_container'
  label 'deepvariant'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/deepvariant"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  tuple path(fa), path(fai), path(dict)
  path bed
  val model_type
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(run), val(dataset), path("*germline.vcf"), emit: germline_vcfs

  script:
  """
  mkdir TMP
  mkdir intermediate_results
  export TMPDIR=\${PWD}/TMP/

  /opt/deepvariant/bin/run_deepvariant \
    --model_type=${model_type} \
    --ref=${fa} \
    --reads=${bam} \
    --regions=\\"${bed}\\" \
    --output_vcf=${dataset}-${pat_name}-${run}${suffix}.germline.vcf \
    --num_shards ${task.cpus} \
    --intermediate_results_dir intermediate_results
  """
}
