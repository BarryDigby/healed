#!/usr/bin/env nextflow

process strelka2_somatic {
// Runs strelka2 on a set of paired samples
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(dataset) - Dataset
//     val(norm_run) - Normal Run Name
//     path(norm_bam) - Normal Alignment File
//     path(norm_bai) - Normal Alignment Index
//     val(tumor_run) - Tumor Run Name
//     path(tumor_bam) - Tumor Alignment File
//     path(tumor_bai) - Tumor Alignment Index
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Reference Index Files
//     path(dict_file) - Reference Dict File
//   parstr - Additional Parameters
//
// output:
//   tuple
//     val(pat_name) - Patient Name
//     val(dataset) - Dataset
//     path('*vcf') - Variant Call File

// require:
//   NORM_TUMOR_BAMS_BAIS
//   REF_W_INDICES
//   params.strelka2$strelka2_somatic_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'strelka2_container'
  label 'strelka2_somatic'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/strelka2_somatic"


  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai)
  tuple path(fa), path(idx_files), path(dict_file)
  tuple path(bed), path(bed_tbi)
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*snvs.vcf.gz"), emit: snv_vcfs
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*indels.vcf.gz"), emit: indel_vcfs

  script:
  """
  /opt/strelka/bin/configureStrelkaSomaticWorkflow.py ${parstr} --tumorBam ${tumor_bam} --normalBam ${norm_bam} --runDir StrelkaSomaticWorkflow --referenceFasta *.fa --exome --callRegions ${bed}
  
  ./StrelkaSomaticWorkflow/runWorkflow.py -j ${task.cpus} -m local
  mv StrelkaSomaticWorkflow/results/variants/* .
  mv somatic.snvs.vcf.gz ${dataset}-${pat_name}-${norm_run}_${tumor_run}${suffix}.snvs.vcf.gz
  mv somatic.indels.vcf.gz ${dataset}-${pat_name}-${norm_run}_${tumor_run}${suffix}.indels.vcf.gz
  """
}
