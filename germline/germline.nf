#!/usr/bin/env nextflow

include { bwa_index } from '../bwa/bwa.nf'
include { bwa_mem } from '../bwa/bwa.nf'
include { gatk_haplotypecaller } from '../gatk4/gatk4.nf'
include { deepvariant } from '../deepvariant/deepvariant.nf'
include { manifest_to_dna_procd_fqs } from '../preproc/preproc.nf'
include {picard_create_seq_dict } from '../picard2/picard2.nf'
include { sams_to_sorted_bams_w_indices } from '../samtools/samtools.nf'
include { samtools_index } from '../samtools/samtools.nf'
include { samtools_faidx } from '../samtools/samtools.nf'
include { trim_galore } from '../trim_galore/trim_galore.nf'
include { combine_sample_files } from '../utilities/utilities.nf'
include { get_fastqs } from '../utilities/utilities.nf'

INDEX_FILES=''

/* HaplotypeCaller workflows under construction */
//workflow manifest_to_haplotypecaller_germline_vcfs {
//// require:
////   MANIFEST
////   params.germline$manifest_to_haplotypecaller_germline_vcfs$dna_reference
////   INDEX_FILES
////   params.germline$manifest_to_haplotypecaller_germline_vcfs$haplotypecaller_parameters
////   params.germline$manifest_to_haplotypecaller_germline_vcfs$haplotypecaller_suffix
//  take:
//    manifest
//    ref
//    idx_files
//    haplotype_caller_parameters
//    suffix
//  main:
//    get_fastqs(
//      manifest,
//      params.fq_dir)
//    raw_fqs_to_haplotypecaller_germline_vcfs(
//      get_fastqs.out.fastqs,
//      ref,
//      idx_files)
//  emit:
//    germline_vcfs = raw_fqs_to_haplotypecaller_germline_vcfs.out.germline_vcfs
//}
//
//
//workflow raw_fqs_to_haplotypecaller_germline_vcfs {
//// require:
////   FQS
////   params.germline$raw_fqs_to_haplotypecaller_germline_vcfs$dna_reference
////   IDX_FILES
//  take:
//    fqs
//    ref
//    idx_files
//  main:
//    trim_galore(
//      fqs,
//      params.germline$raw_fqs_to_haplotypecaller_germline_vcfs$trim_galore_parameters)
//    procd_fqs_to_haplotypecaller_germline_vcfs(
//      trim_galore.out.procd_fqs,
//      ref,
//      idx_files)
//  emit:
//    germline_vcfs = procd_fqs_to_haplotypecaller_germline_vcfs.out.germline_vcfs
//}
//
//
//workflow procd_fqs_to_haplotypecaller_germline_vcfs {
//  take:
//    fqs
//    ref
//    idx_files
//  main:
//    if( idx_files == '' ) {
//      bwa_index(
//        ref,
//        params.germline$procd_fqs_to_haplotypecaller_germline_vcfs$bwa_index_parameters)
//        bwa_index.out.idx_files.set{ idx_files }
//    }
//    bwa_mem_samtools_sort(
//      fqs,
//      idx_files,
//      params.germline$procd_fqs_to_haplotypecaller_germline_vcfs$bwa_mem_parameters)
//    raw_alns_to_haplotypecaller_germline_vcfs(
//      bwa_mem_samtools_sort.out.bams,
//      ref)
//  emit:
//    germline_vcfs = raw_alns_to_haplotypecaller_germline_vcfs.out.germline_vcfs
//}
//
//
//workflow raw_alns_to_haplotypecaller_germline_vcfs {
//  take:
//    bams
//    ref
//  main:
//    sams_to_sorted_bams_w_indices(
//      bams)
//    picard_mark_duplicates(
//      sams_to_sorted_bams_w_indices,
//      params.germline$raw_alns_to_haplotypecaller_germline_vcfs$picard_mark_duplicates_parameters)
//    samtools_index(
//      picard_mark_duplicates.out.mkdup_bams,
//      '')
//    combine_sample_files(picard_mark_duplicates.out.mkdup_bams, samtools_index.out.bais).set{ bams_and_bais }
//    procd_alns_to_haplotypecaller_germline_vcfs(
//      bams_and_bais,
//      ref)
//  emit:
//    germline_vcfs = procd_alns_to_haplotypecaller_germline_vcfs.out.germline_vcfs
//}
//
//
//workflow procd_alns_to_haplotypecaller_germline_vcfs {
//  take:
//    bams
//    ref
//  main:
//    make_ancillary_index_files(
//      ref)
//    gatk_haplotypecaller(
//      bams,
//      make_ancillary_index_files.out.collective_idx_files,
//      params.germline$procd_alns_to_haplotype_germline_vcfs$haplotypecaller_parameters,
//      params.germline$procd_alns_to_haplotype_germline_vcfs$haplotypecaller_suffix)
//  emit:
//    germline_vcfs = gatk_haplotypecaller.out.germline_vcfs
//}


//Maybe this shouldn't be in somatic, but rather in a generic alignments tools module?
//I think this is identical is get_ref_idx_files
workflow make_ancillary_index_files {
// take:
//   ref - Reference FASTA
//   samtools_faidx_parameters - SAMtools faidx Parameters
//
// emit:

// require:
//   params.germline$make_ancillary_index_files$dna_Ref
//   params.germline$make_ancillary_index_files$samtools_faidx_parameters
  take:
    ref
    samtools_faidx_parameters
  main:
    samtools_faidx(
      ref,
      samtools_faidx_parameters)
    picard_create_seq_dict(
      ref,
      '')
    //Can't join since the ref symlinks are different.
    samtools_faidx.out.faidx_file
      .concat(picard_create_seq_dict.out.dict_file)
      .collect().map{ [it[0], it[1], it[3]] }
      .set{ collective_idx_files }
  emit:
    collective_idx_files
}


workflow manifest_to_deepvariant_germline_vcfs {
// require:
//   MANIFEST
//   IDX_FILES
//   params.germline$manifest_to_deepvariant_germline_vcfs$dna_ref
//   params.germline$manifest_to_deepvariant_germline_vcfs$trim_galore_parameters
//   params.germline$manifest_to_deepvariant_germline_vcfs$bwa_index_parameters
//   params.germline$manifest_to_deepvariant_germline_vcfs$bwa_mem_norm_parameters
//   params.germline$manifest_to_deepvariant_germline_vcfs$picard_mark_duplicates_parameters
//   params.germline$manifest_to_deepvariant_germline_vcfs$samtools_index_parameters
//   params.germline$manifest_to_deepvariant_germline_vcfs$deepvariant_model_type
//   params.germline$manifest_to_deepvariant_germline_vcfs$deepvariant_parameters
//   params.germline$manifest_to_deepvariant_germline_vcfs$deepvariant_suffix
  take:
    manifest
    idx_files
    dna_ref
    trim_galore_parameters 
    bwa_index_parameters
    bwa_mem_norm_parameters
    picard_mark_duplicates_parameters
    samtools_index_parameters
    deepvariant_model_type
    deepvariant_parameters
    deepvariant_suffix
  main:
    get_fastqs(
      manifest,
      params.fq_dir)
    raw_fqs_to_deepvariant_germline_vcfs(
      get_fastqs.out.fastqs,
      idx_files,
      dna_ref,
      trim_galore_parameters,
      bwa_index_parameters,
      bwa_mem_norm_paramters,
      picard_mark_duplicates_parameters,
      samtools_index_parameters,
      targets_bed,
      deepvariant_model_type,
      deepvariant_parameters,
      deepvariant_suffix)
  emit:
    germline_vcfs = raw_fqs_to_deepvariant_germline_vcfs.out.germline_vcfs
}


workflow raw_fqs_to_deepvariant_germline_vcfs {
// require:
//   FQS
//   IDX_FILES
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$dna_ref
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$trim_galore_parameters
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$bwa_index_parameters
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$bwa_mem_norm_parameters
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$picard_mark_duplicates_parameters
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$samtool_index_parameters
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$targets_bed
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$deepvariant_model_type
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$deepvariant_parameters
//   params.germline$raw_fqs_to_deepvariant_germline_vcfs$deepvariant_suffix
  take:
    fqs
    idx_files
    dna_ref
    trim_galore_parameters 
    bwa_index_parameters
    bwa_mem_norm_parameters
    picard_mark_duplicates_parameters
    samtools_index_parameters
    targets_bed
    deepvariant_model_type
    deepvariant_parameters
    deepvariant_suffix
  main:
    trim_galore(
      fqs,
      trim_galore_parameters)
    procd_fqs_to_deepvariant_germline_vcfs(
      trim_galore.out.procd_fqs,
      idx_files,
      dna_ref,
      bwa_index_parameters,
      bwa_mem_norm_parameters,
      picard_mark_duplicates_parameters,
      samtools_index_parameters,
      targets_bed,
      deepvariant_model_type,
      deepvariant_parameters,
      deepvariant_suffix)
  emit:
    germline_vcfs = procd_fqs_to_deepvariant_germline_vcfs.out.germline_vcfs
}


workflow procd_fqs_to_deepvariant_germline_vcfs {
// require
//   FQS
//   IDX_FILES
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$dna_ref
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$bwa_index_parameters
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$bwa_mem_norm_parameters
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$picard_mark_duplicates_parameters
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$samtools_index_parameters
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$targets_bed
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$deepvariant_model_type
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$deepvariant_parameters
//   params.germline$procd_fqs_to_deepvariant_germline_vcfs$deepvariant_suffix
  take:
    fqs
    idx_files
    dna_ref
    bwa_index_parameters
    bwa_mem_norm_parameters
    picard_mark_duplicates_parameters
    samtools_index_parameters
    targets_bed
    deepvariant_model_type
    deepvariant_parameters
    deepvariant_suffix
  main:
    if( idx_files == '' ) {
      bwa_index(
        dna_ref,
        bwa_index_parameters)
        bwa_index.out.idx_files.set{ idx_files }
    }
    bwa_mem_samtools_sort(
      fqs,
      idx_files,
      bwa_mem_norm_parameters)
    raw_alns_to_deepvariant_germline_vcfs(
      bwa_mem_samtools_sort.out.bams,
      dna_ref,
      picard_mark_duplicates_parameters,
      samtools_index_parameters,
      targets_bed,
      deepvariant_model_type,
      deepvariant_parameters,
      deepvariant_suffix)
  emit:
    germline_vcfs = raw_alns_to_deepvariant_germline_vcfs.out.germline_vcfs
}


workflow raw_alns_to_deepvariant_germline_vcfs {
// require:
//   BAMS
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$dna_ref
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$picard_mark_duplicates_parameters
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$samtools_index_parameters
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$targets_bed
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$deepvariant_model_type
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$deepvariant_parameters
//   params.germline$raw_alns_to_deepvariant_germline_vcfs$deepvariant_suffix
  take:
    bams
    dna_ref
    picard_mark_duplicates_parameters
    samtools_index_parameters
    targets_bed
    deepvariant_model_type
    deepvariant_parameters
    deepvariant_suffix
  main:
    sams_to_sorted_bams_w_indices(
      bams)
    picard_mark_duplicates(
      sams_to_sorted_bams_w_indices,
      picard_mark_duplicates_parameters)
    samtools_index(
      picard_mark_duplicates.out.mkdup_bams,
      samtools_index_parameters)
    combine_sample_files(
      picard_mark_duplicates.out.mkdup_bams,
      samtools_index.out.bais)
      .set{ bams_and_bais }
    procd_alns_to_deepvariant_germline_vcfs(
      bams_and_bais,
      dna_ref,
      targets_bed,
      deepvariant_model_type,
      deepvariant_parameters,
      deepvariant_suffix)
  emit:
    germline_vcfs = procd_alns_to_deepvariant_germline_vcfs.out.germline_vcfs
}


workflow procd_alns_to_deepvariant_germline_vcfs {
// require:
//   BAMS
//   params.germline$procd_alns_to_deepvariant_germline_vcfs$dna_ref
//   params.germline$procd_alns_to_deepvariant_germline_vcfs$targets_bed
//   params.germline$procd_alns_to_deepvariant_germline_vcfs$deepvariant_model_type
//   params.germline$procd_alns_to_deepvariant_germline_vcfs$deepvariant_parameters
//   params.germline$procd_alns_to_deepvariant_germline_vcfs$deepvariant_suffix
  take:
    bams
    dna_ref
    targets_bed
    deepvariant_model_type
    deepvariant_parameters
    deepvariant_suffix
    samtools_faidx_parameters
  main:
    make_ancillary_index_files(
      dna_ref,
      samtools_faidx_parameters)
    deepvariant(
      bams,
      make_ancillary_index_files.out.collective_idx_files,
      targets_bed,
      deepvariant_model_type,
      deepvariant_parameters,
      deepvariant_suffix)
  emit:
    germline_vcfs = deepvariant.out.germline_vcfs
}
