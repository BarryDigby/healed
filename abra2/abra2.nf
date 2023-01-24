#!/usr/bin/env nextflow

include { combine_sample_files }    from '../utilities/utilities.nf'
include { combine_patient_samples } from '../utilities/utilities.nf'


process abra2_cadabra {
// Runs abra2.jar abra.cadabra.Cadabra
// Cadabra is a somatic indel caller that works specifically with ABRA alignments.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(dataset) - Dataset
//     val(norm_run) - Normal Run Name
//     path(norm_bam) - Normal BAM
//     path(norm_bai) - Normal BAI
//     val(tumor_run) - Tumor Run Name
//     path(tumor_bam) - Tumor BAM
//     path(tumor_bai) - Tumor BAI
//   path fa - Reference FASTA
//   val suffix - Output VCF Suffix
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: vcfs
//     val(pat_name) -  Patient Name
//     val(norm_run) -  Normal Run Name
//     val(tumor_run) -  Tumor Run Name
//     val(dataset) -  Dataset
//     path("*.vcf") - Output VCF File

// require:
//   NORM_TUM_BAMS_AND_BAIS
//   REF
//   params.abra2$abra2_cadabra_suffix
//   params.abra2$abra2_cadabra_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'abra2_container'
  label 'abra2_cadabra'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/abra2_cadabra"

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai)
  path fa
  val suffix
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.vcf"), emit: vcfs

  script:
  """
  java -Xmx${task.memory.toGiga()}g -cp /abra2.jar abra.cadabra.Cadabra --threads ${task.cpus} --ref ${fa} --normal ${norm_bam} --tumor ${tumor_bam} ${parstr} > ${dataset}-${pat_name}-${norm_run}_${tumor_run}${suffix}.vcf
  """
}


process abra2 {
// Runs abra2
// ABRA is a realigner for next generation sequencing data. It uses localized
// assembly and global realignment to align reads more accurately, thus
// improving downstream analysis (detection of indels and complex variants in
// particular).
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(datset) - Dataset
//     val(norm_run) - Normal Run Name
//     val(norm_bam) - Normal BAM
//     val(norm_bai) - Normal BAI
//     val(tumor_run) - Tumor Run Name
//     val(tumor_bam) - Tumor BAM
//     path(tumor_bai) - Tumor BAI
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Reference Index Files
//   val bed - BED File
//   val parstr - Additional Parameters
//   val tmp_dir - Temporary Directory
//
// output:
//   tuple => emit: norm_abra_bams
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(dataset) - Dataset
//     path("${norm_run}*abra.bam") - Normal Realigned Output BAM
//   tuple => emit: tumor_abra_bams
//     val(pat_name) - Patient Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("${tumor_run}*abra.bam") - Tumor Realigned Output BAM
//   tuple => emit: norm_tumor_abra_bams
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("${norm_run}*abra.bam") - Normal Realigned Output BAM
//     path("${tumor_run}*abra.bam") - Tumor Realigned Output BAM

// require:
//   ALNS
//   REF_W_INDEX
//   params.abra2$abra2_parameters
//   params.abra2$abra2_tmp_dir

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'abra2_container'
  label 'abra2_realign'

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai)
  tuple path(fa), path(idx_files)
  path bed
  val parstr
  val tmp_dir

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*norm_abra.bam"), emit: norm_abra_bams
  tuple val(pat_name), val(tumor_run), val(norm_run), val(dataset), path("*tumor_abra.bam"), emit: tumor_abra_bams
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*norm_abra.bam"), path("*tumor_abra.bam"), emit: norm_tumor_abra_bams

  script:
  """
  out_norm_bam_bfr=`echo ${norm_bam}`
  out_norm_bam=\${out_norm_bam_bfr%.bam}.norm_abra.bam

  out_tumor_bam_bfr=`echo ${tumor_bam}`
  out_tumor_bam=\${out_tumor_bam_bfr%.bam}.tumor_abra.bam

  mkdir -p ${tmp_dir}

  java -Xmx${task.memory.toGiga()}g -jar /abra2.jar --in ${norm_bam},${tumor_bam} --out \${out_norm_bam},\${out_tumor_bam} --ref ${fa} --targets ${bed} --tmpdir ${tmp_dir} --threads ${task.cpus} ${parstr} > abra.log
  """
}


process abra2_rna {
// Runs abra2 on RNA BAMs
// ABRA is a realigner for next generation sequencing data. It uses localized
// assembly and global realignment to align reads more accurately, thus
// improving downstream analysis (detection of indels and complex variants in
// particular).
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(datset) - Dataset
//     val(run) - RNA Run Name
//     val(norm_bam) - Normal BAM
//     val(norm_bai) - Normal BAI
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Reference Index Files
//   path gtf - Reference GTF File
//   val parstr - Additional Parameters
//   val tmp_dir - Temporary Directory
//
// output:
//   tuple => emit: abra_bams
//     val(pat_name) - Patient Name
//     val(run) - RNA Run Name
//     val(dataset) - Dataset
//     path("${run}*abra.bam") - RNA Realigned Output BAM

// require:
//   ALNS
//   REF_W_INDEX
//   params.abra2$abra2_parameters
//   params.abra2$abra2_tmp_dir

  tag "${dataset}/${pat_name}/${run}"
  label 'abra2_container'
  label 'abra2_realign_rna'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path(junctions), path(vcf)
  tuple path(fa), path(idx_files)
  path gtf
  path bed
  val parstr
  val tmp_dir

  output:
  tuple val(pat_name), val(run), val(dataset), path("*${run}*abra.bam"), emit: abra_bams

  script:
  """
  bam_prec=`echo ${bam}`
  out_bam=\${bam_prec%.bam}.abra.bam

  mkdir ${tmp_dir}

  java -Xmx${task.memory.toGiga()}g -jar /abra2.jar --in ${bam} --junctions ${junctions} --in-vcf ${vcf} --out \${out_bam} --ref ${fa} --targets ${bed} --gtf ${gtf} --threads ${task.cpus} --tmpdir ${tmp_dir} ${parstr} > abra.log
  """
}


workflow bams_bais_to_realigned_bams {
// Produces ABRA2 realigned BAMS given reference FASTA, BAMs, BAIs, and
// manifest channels
//
// NOTE: There should be a 1-to-1-to-1 mapping among samples, BAMs, and BAIs
// among the manifest, bams, and bais channels.
//
// take:
//   dna_ref - Reference FASTA
//   targets_bed - Targets BED file (e.g. exome)
//   abra2_parameters - Abra2 Parameters
//   bams - BAMs
//   bais - BAIs
//   manifest - RAFT Manifest
//
// emit:
//    norm_abra_bams - Normal Realigned BAMs
//    tumor_abra_bams - Tumor Realigned BAMS
//    norm_tumor_abra_bams - Both Normal and Tumor Realigned BAMs

// require:
//   params.abra2$bams_bais_to_realigned_bams$dna_ref
//   params.abra2$bams_bais_to_realigned_bams$targets_bed,
//   params.abra2$bams_bais_to_realigned_bams$abra2_parameters,
//   BAMS
//   BAIS
//   MANIFEST

  take:
    dna_ref
    targets_bed
    abra2_parameters
    bams
    bais
    manifest
  main:
    combine_sample_files(
      bams,
      bais)
    combine_patient_samples(
      combine_sample_files.out.combined_set,
      'TRUE', //Normal sample
      'FALSE', //Abnormal sample
      manifest)
    abra2(
      combine_patient_samples.out.combined_set,
      dna_ref,
      targets_bed,
      abra2_parameters,
      'tmp_dir')
    // This will be explained further in the wiki, but wanted to include some
    // context to the "remapped" channels. Most tools that take a "normal" and
    // a "tumor" input combine them in such a way to make a single resulting
    // output (e.g. a VCF). ABRA is unique in that it creates two separate
    // outputs, a normal BAM realigned using the tumor bam and a tumor BAM
    // realigned using the normal bam. One could use the paired run (e.g.
    // val(norm_run) and val(tumor_run) to keep these outputs linked, but
    // many of the downstream processses and workflows expect to have a single
    // run. As a result, the sample-level run is remapped to contain both
    // normal and tumor runs information. This is later deconvulated
    // downstream.
    abra2.out.norm_abra_bams.map{ [it[0], "${it[1]}-rel-${it[1]}-to-${it[2]}", it[3], it[4]] }.set{ remapped_norm_abra_bams }
    abra2.out.tumor_abra_bams.map{ [it[0], "${it[1]}-rel-${it[2]}-to-${it[1]}", it[3], it[4]] }.set{ remapped_tumor_abra_bams }
    abra2.out.norm_tumor_abra_bams
      .flatMap{ [[it[0], "${it[1]}-rel-${it[1]}-to-${it[2]}", it[3], it[4]], [it[0], "${it[2]}-rel-${it[1]}-to-${it[2]}", it[3], it[5]]] }
      .set{ remapped_abra_bams }
    remapped_abra_bams.set{ norm_tumor_abra_bams }
  emit:
    norm_tumor_abra_bams = norm_tumor_abra_bams
    norm_abra_bams = abra2.out.norm_abra_bams
    tumor_abra_bams = abra2.out.tumor_abra_bams
}
