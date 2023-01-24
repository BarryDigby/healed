#!/usr/bin/env nextflow


include { abra2_cadabra } from '../abra2/abra2.nf'
include { bams_bais_to_realigned_bams } from '../abra2/abra2.nf'

include { bedtools_intersect } from '../bedtools/bedtools.nf'

include { bwa_index } from '../bwa/bwa.nf'
include { bwa_mem_samtools_sort as bwa_mem } from '../bwa/bwa.nf'

include { gatk_apply_bqsr } from '../gatk4/gatk4.nf'
include { gatk_base_recalibrator } from '../gatk4/gatk4.nf'
include { gatk_index_feature_file } from '../gatk4/gatk4.nf'
include { gatk_mutect2_matched } from '../gatk4/gatk4.nf'

include { picard_create_seq_dict } from '../picard2/picard2.nf'
include { picard_collect_insert_size_metrics } from '../picard2/picard2.nf'
include { picard_mark_duplicates } from '../picard2/picard2.nf'

include { manifest_to_dna_procd_fqs } from '../preproc/preproc.nf'

include { bams_to_sorted_bams_w_indices } from '../samtools/samtools.nf'
include { samtools_faidx } from '../samtools/samtools.nf'
include { samtools_index } from '../samtools/samtools.nf'
include { samtools_index as samtools_index_recal } from '../samtools/samtools.nf'
include { samtools_index as samtools_index_bqsr } from '../samtools/samtools.nf'
include { samtools_rmdup } from '../samtools/samtools.nf'
include { samtools_sort } from '../samtools/samtools.nf'
include { samtools_view } from '../samtools/samtools.nf'

include { htslib_bgzip } from '../htslib/htslib.nf'
include { htslib_tabix } from '../htslib/htslib.nf'
include { htslib_tabix as htslib_tabix_pon } from '../htslib/htslib.nf'
include { htslib_tabix as htslib_tabix_af } from '../htslib/htslib.nf'

include { get_callable_bed } from '../ref_utils/ref_utils.nf'

include { strelka2_somatic } from '../strelka2/strelka2.nf'

include { trim_galore } from '../trim_galore/trim_galore.nf'

include { get_fastqs } from '../utilities/utilities.nf'
include { filter_channel_by_manifest as filter_channel_by_manifest_norm } from '../utilities/utilities.nf'
include { filter_channel_by_manifest as filter_channel_by_manifest_tumor } from '../utilities/utilities.nf'
include { combine_patient_samples } from '../utilities/utilities.nf'

//Defaulting to null string to prompt creation of index files if needed.
INDEX_FILES=''


workflow manifest_to_mutect2_somatic_vcfs {
// require:
//   params.somatic$manifest_to_mutect2_somatic_vcfs$dna_ref
//   INDEX_FILES
//   MANIFEST
  take:
    ref
    index_files
    manifest
  main:
    get_fastqs(
      manifest,
      params.somatic$fq_dir)
    raw_fqs_to_mutect2_somatic_vcfs(
      ref,
      index_files,
      get_fastqs.out.fastqs,
      manifest)
  emit:
    fqs = get_fastqs.out.fastqs
    procd_fqs = raw_fqs_to_mutect2_somatic_vcfs.out.procd_fqs
    raw_alns = raw_fqs_to_mutect2_somatic_vcfs.out.raw_alns
    norm_tumor_procd_bams_bais = raw_fqs_to_mutect2_somatic_vcfs.out.norm_tumor_procd_bams_bais
    vcfs = raw_fqs_to_mutect2_somatic_vcfs.out.vcfs
}


workflow manifest_to_cadabra_somatic_vcfs {
// require:
//   params.somatic$manifest_to_cadabra_somatic_vcfs$dna_ref
//   INDEX_FILES
//   MANIFEST
  take:
    ref
    index_files
    manifest
  main:
    get_fastqs(
      manifest,
      params.somatic$fq_dir)
    raw_fqs_to_cadabra_somatic_vcfs(
      ref,
      index_files,
      get_fastqs.out.fastqs,
      manifest)
  emit:
    vcfs = raw_fqs_to_cadabra_somatic_vcfs.out.vcfs
}


workflow manifest_to_strelka2_somatic_vcfs {
// require:
//   params.somatic$manifest_to_strelka2_somatic_vcfs$dna_ref
//   INDEX_FILES
//   MANIFEST
  take:
    ref
    index_files
    manifest
  main:
    get_fastqs(
      manifest,
      params.somatic$fq_dir)
    raw_fqs_to_strelka2_somatic_vcfs(
      ref,
      index_files,
      get_fastqs.out.fastqs,
      manifest)
  emit:
    vcfs = raw_fqs_to_strelka2_somatic_vcfs.out.vcfs
}


workflow raw_fqs_to_mutect2_somatic_vcfs {
// require:
//   params.somatic$raw_fqs_to_mutect2_somatic_vcfs$dna_ref
//   INDEX_FILES
//   FQS
//   MANIFEST
  take:
    ref
    index_files
    fqs
    manifest
  main:
    trim_galore(
      fqs,
      params.somatic$raw_fqs_to_mutect2_somatic_vcfs$trim_galore_parameters)
    procd_fqs_to_mutect2_somatic_vcfs(
      ref,
      index_files,
      trim_galore.out.procd_fqs,
      manifest)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    raw_alns = procd_fqs_to_mutect2_somatic_vcfs.out.raw_alns
    norm_tumor_procd_bams_bais = procd_fqs_to_mutect2_somatic_vcfs.out.norm_tumor_procd_bams_bais
    vcfs = procd_fqs_to_mutect2_somatic_vcfs.out.vcfs
}


workflow raw_fqs_to_cadabra_somatic_vcfs {
// require:
//   params.somatic$raw_fqs_to_cadabra_somatic_vcfs$dna_ref
//   INDEX_FILES
//   FQS
//   MANIFEST
  take:
    ref
    index_files
    fqs
    manifest
  main:
    trim_galore(
      fqs,
      params.somatic$raw_fqs_to_cadabra_somatic_vcfs$trim_galore_parameters)
    procd_fqs_to_cadabra_somatic_vcfs(
      ref,
      index_files,
      trim_galore.out.procd_fqs,
      manifest)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    raw_alns = procd_fqs_to_cadabra_somatic_vcfs.out.raw_alns
    norm_tumor_procd_bams_bais = procd_fqs_to_cadabra_somatic_vcfs.out.norm_tumor_procd_bams_bais
    vcfs = procd_fqs_to_cadabra_somatic_vcfs.out.vcfs
}


workflow raw_fqs_to_strelka2_somatic_vcfs {
// require:
//   params.somatic$raw_fqs_to_strelka2_somatic_vcfs$dna_ref
//   INDEX_FILES
//   FQS
//   MANIFEST
  take:
    ref
    index_files
    fqs
    manifest
  main:
    trim_galore(
      fqs,
      params.somatic$raw_fqs_to_strelka2_somatic_vcfs$trim_galore_parameters)
    procd_fqs_to_strelka2_somatic_vcfs(
      ref,
      index_files,
      trim_galore.out.procd_fqs,
      manifest)
  emit:
    vcfs = procd_fqs_to_strelka2_somatic_vcfs.out.vcfs
}


workflow procd_fqs_to_mutect2_somatic_vcfs {
// require:
//   params.somatic$procd_fqs_to_mutect2_somatic_vcfs$dna_ref
//   INDEX_FILES
//   FQS
//   MANIFEST
  take:
    ref
    index_files
    fqs
    manifest
  main:
    if( index_files == '' ) {
      bwa_index(
        ref,
        params.somatic$procd_fqs_to_mutect2_somatic_vcfs$bwa_index_parameters)
      bwa.index.out.idx_files.set{ index_files }
    }
    bwa_mem(
      fqs,
      index_files,
      params.somatic$procd_fqs_to_mutect2_somatic_vcfs$bwa_mem_parameters)
    raw_alns_to_mutect2_somatic_vcfs(
      ref,
      index_files,
      bwa_mem.out.sams,
      manifest)
  emit:
    raw_alns = bwa_mem.out.sams
    norm_tumor_procd_bams_bais = raw_alns_to_mutect2_somatic_vcfs.out.norm_tumor_procd_bams_bais
    vcfs = raw_alns_to_mutect2_somatic_vcfs.out.vcfs
}


workflow procd_fqs_to_cadabra_somatic_vcfs {
// require:
//   params.somatic$procd_fqs_to_cadabra_somatic_vcfs$dna_ref
//   INDEX_FILES
//   FQS
//   MANIFEST
  take:
    ref
    index_files
    fqs
    manifest
  main:
    if( index_files == '' ) {
      bwa_index(ref,
        params.somatic$procd_fqs_to_cadabra_somatic_vcfs$bwa_index_parameters)
      bwa.index.out.idx_files.set{ index_files }
    }
    bwa_mem(
      fqs,
      index_files,
      params.somatic$procd_fqs_to_cadabra_somatic_vcfs$bwa_mem_parameters)
    raw_alns_to_cadabra_somatic_vcfs(
      ref,
      index_files,
      bwa_mem.out.sams,
      manifest)
  emit:
    vcfs = raw_alns_to_cadabra_somatic_vcfs.out.vcfs
}


workflow procd_fqs_to_strelka2_somatic_vcfs {
// require:
//   params.somatic$procd_fqs_to_strelka2_somatic_vcfs$dna_ref
//   INDEX_FILES
//   FQS
//   MANIFEST
  take:
    ref
    index_files
    fqs
    manifest
  main:
    if( index_files == '' ) {
      bwa_index(
        dna_ref,
        params.somatic$procd_fqs_to_strelka2_somatic_vcfs$bwa_index_parameters)
      bwa.index.out.idx_files.set{ index_files }
    }
    bwa_mem(
      fqs,
      index_files,
      params.somatic$procd_fqs_to_strelka2_somatic_vcfs$bwa_mem_parameters)
    raw_alns_to_strelka2_somatic_vcfs(
      ref,
      index_files,
      bwa_mem.out.sams,
      manifest)
  emit:
    vcfs = raw_alns_to_strelka2_somatic_vcfs.out.vcfs
}


//Should this be here or in preproc.nf?
workflow procd_fqs_to_raw_alns {
// require:
//   FQS
//   MANIFEST
//   params.somatic$procd_fqs_to_raw_alns$dna_ref
//   params.somatic$procd_fqs_to_raw_alns$bwa_index_parameters
//   params.somatic$procd_fqs_to_raw_alns$bwa_mem_parameters
  take:
    fqs
    manifest
    dna_ref
    bwa_index_parameters
    bwa_mem_parameters
  main:
    bwa_index(
      dna_ref,
      bwa_index_parameters)
    bwa_mem(
      fqs,
      bwa_index.out.idx_files,
      bwa_mem_parameters)
  emit:
    alns = bwa_mem.out.bams
    idx_files = bwa_index.out.idx_files
}


//Should this be here or in preproc.nf?
workflow procd_fqs_to_procd_alns {
// require:
//    FQS
//    MANIFEST
//    params.somatic$procd_fqs_to_procd_alns$dna_ref
  take:
    fqs
    manifest
    dna_ref
    bwa_index_parameters
    bwa_mem_parameters
    targets_bed
    picard_mark_duplicates_parameters
    known_sites_ref
    samtools_index_parameters
    gatk_index_feature_file_parameters
    gatk_base_recalibrator_parameters
    gatk_apply_bqsr_parameters
    samtools_faidx_parameters
    abra2_parameters
  main:
    procd_fqs_to_raw_alns(
      fqs,
      manifest,
      dna_ref,
      bwa_index_parameters,
      bwa_mem_parameters)
    raw_alns_to_procd_alns(
      dna_ref,
      targets_bed,
      picard_mark_duplicates_parameters,
      known_sites_ref,
      procd_fqs_to_raw_alns.out.alns,
      procd_fqs_to_raw_alns.out.idx_files,
      manifest,
      samtools_index_parameters,
      gatk_index_feature_file_parameters,
      gatk_base_recalibrator_parameters,
      gatk_apply_bqsr_parameters,
      samtools_faidx_parameters,
      abra2_parameters)
  emit:
    filt_bams_and_bais = raw_alns_to_procd_alns.out.filt_bams_and_bais
    filt_bams = raw_alns_to_procd_alns.out.filt_bams
    filt_bais = raw_alns_to_procd_alns.out.filt_bais
    anc_idx_files = raw_alns_to_procd_alns.out.anc_idx_files
    marked_dup_metrics = raw_alns_to_procd_alns.out.marked_dup_metrics
}


workflow raw_alns_to_mutect2_somatic_vcfs {
// require:
//   params.somatic$raw_alns_to_mutect2_somatic_vcfs$dna_ref
//   params.somatic$raw_alns_to_mutect2_somatic_vcfs$targets_bed
//   params.somatic$raw_alns_to_mutect2_somatic_vcfs$abra2_parameters
//   params.somatic$raw_alns_to_mutect2_somatic_vcfs$idx_files
//   params.somatic$raw_alns_to_mutect2_somatic_vcfs$alns
//   params.somatic$raw_alns_to_mutect2_somatic_vcfs$manifest
  take:
    dna_ref
    targets_bed
    abra2_parameters
    idx_files
    alns
    manifest
  main:
    samtools_faidx(
      dna_ref,
      params.somatic$indices_out_dir,
      params.shared_dir)
    picard_create_seq_dict(
      dna_ref,
      params.somatic$indices_out_dir,
      params.shared_dir)
    //Can't join since the ref symlinks are different.
    samtools_faidx.out.faidx_file
      .concat(picard_create_seq_dict.out.dict_file)
      .collect()
      .map{ [it[0], it[1], it[3]] }
      .set{ collective_idx_files }
    samtools_view(
      alns,
      params.somatic$raw_alns_to_mutect2_vcf$samtools_view_params,
      params.somatic$samps_out_dir,
      params.shared_dir)
    samtools_sort(
      samtools_view.out.bams,
      params.somatic$samps_out_dir,
      params.shared_dir)
    samtools_index(
      samtools_sort.out.bams,
      '')
    bams_bais_to_realigned_bams(
      dna_ref,
      targets_bed,
      abra2_parameters,
      samtools_sort.out.bams,
      samtools_index.out.bais,
      manifest)
    picard_collect_insert_size_metrics(
      bams_bais_to_realigned_bams.out.norm_tumor_abra_bams,
      params.somatic$samps_out_dir,
      params.shared_dir)
    bams_bais_to_realigned_bams.out.norm_tumor_abra_bams
      .join(picard_collect_insert_size_metrics.out.insert_size_metrics, by: [0, 1, 2])
      .set{ norm_tumor_abra_bams_w_is_metrics }
    index_feature_file(
      params.somatic$known_sites,
      params.somatic$indices_out_dir,
      params.shared_dir)
  emit:
    norm_tumor_procd_bams_bais = norm_tumor_procd_bams_bais
    vcfs = procd_alns_to_mutect2_somatic_vcfs.out.vcfs
}


workflow get_ref_idx_set {
  take:
    ref
  main:
    samtools_faidx(
      ref,
      params.somatic$indices_out_dir,
      params.shared_dir)
    picard_create_seq_dict(
      ref,
      params.somatic$indices_out_dir,
      params.shared_dir)
    //Can't join since the ref symlinks are different.
    samtools_faidx.out.faidx_files
      .concat(picard_create_seq_dict.out.dict_files)
      .collect()
      .map{ [it[0], it[1], it[3]] }
      .set{ idx_files_set }
  emit:
    idx_files_set
}


workflow raw_alns_to_procd_alns {
//Mimics GATK best practices. Can probably be improved.
// require:
//   params.somatic$raw_alns_to_procd_alns$dna_ref
//   params.somatic$raw_alns_to_procd_alns$targets_bed
//   params.somatic$raw_alns_to_procd_alns$picard_mark_duplicates_parameters
//   params.somatic$raw_alns_to_procd_alns$known_sites_ref
//   ALNS
//   IDX_FILES
//   MANIFEST
  take:
    dna_ref
    targets_bed
    picard_mark_duplicates_parameters
    known_sites_ref
    alns
    idx_files
    manifest
    samtools_index_parameters
    gatk_index_feature_file_parameters
    gatk_base_recalibrator_parameters
    gatk_apply_bqsr_parameters
    samtools_faidx_parameters
    abra2_parameters
  main:
    make_ancillary_index_files(
      dna_ref,
      samtools_faidx_parameters)
    bams_to_sorted_bams_w_indices(
      alns,
      samtools_index_parameters)
    bams_bais_to_realigned_bams(
      idx_files,
      targets_bed,
      abra2_parameters,
      bams_to_sorted_bams_w_indices.out.bams,
      bams_to_sorted_bams_w_indices.out.bais,
      manifest)
    picard_mark_duplicates(
      bams_bais_to_realigned_bams.out.norm_tumor_abra_bams,
      picard_mark_duplicates_parameters)
    bams_to_base_qual_recal_w_indices(
      picard_mark_duplicates.out.mkdup_bams,
      known_sites_ref,
      make_ancillary_index_files.out.collective_idx_files,
      samtools_index_parameters,
      gatk_index_feature_file_parameters,
      targets_bed,
      gatk_base_recalibrator_parameters,
      gatk_apply_bqsr_parameters)
  emit:
    filt_bams_and_bais = bams_to_base_qual_recal_w_indices.out.bams_and_bais
    filt_bams = bams_to_base_qual_recal_w_indices.out.bams
    filt_bais = bams_to_base_qual_recal_w_indices.out.bais
    anc_idx_files = make_ancillary_index_files.out.collective_idx_files
    marked_dup_metrics = picard_mark_duplicates.out.marked_dup_metrics
}


workflow raw_alns_to_cadabra_somatic_vcfs {
// require:
//   params.somatic$raw_alns_to_cadabra_somatic_vcfs$dna_ref
//   params.somatic$raw_alns_to_cadabra_somatic_vcfs$targets_bed
//   params.somatic$raw_alns_to_cadabra_somatic_vcfs$abra2_parameters
//   params.somatic$raw_alns_to_cadabra_somatic_vcfs$idx_files
//   params.somatic$raw_alns_to_cadabra_somatic_vcfs$alns
//   params.somatic$raw_alns_to_cadabra_somatic_vcfs$manifest
  take:
    dna_ref
    targets_bed
    abra2_parameters
    idx_files
    alns
    manifest
  main:
    samtools_faidx(
      dna_ref,
      params.somatic$indices_out_dir,
      params.shared_dir)
    picard_create_seq_dict(
      dna_ref,
      params.somatic$indices_out_dir,
      params.shared_dir)
    //Can't join since the ref symlinks are different.
    samtools_faidx.out.faidx_file
      .concat(picard_create_seq_dict.out.dict_file)
      .collect()
      .map{ [it[0], it[1], it[3]] }
      .set{ collective_idx_files }
    sams_to_sorted_bams_w_indices(
      alns,
      params.somatic$raw_alns_to_cadabra_somatic_vcfs$samtools_view_params,
      params.somatic$samps_out_dir,
      params.shared_dir)
    bams_bais_to_realigned_bams(
      idx_files,
      targets_bed,
      abra2_parameters,
      sams_to_sorted_bams_w_indices.out.bams,
      sams_to_sorted_bams_w_indices.out.bais,
      manifest)
    picard_collect_insert_size_metrics(
      bams_bais_to_realigned_bams.out.norm_tumor_abra_bams,
      params.somatic$samps_out_dir,
      params.shared_dir)
    bams_bais_to_realigned_bams.out.norm_tumor_abra_bams
      .join(picard_collect_insert_size_metrics.out.insert_size_metrics, by: [0, 1, 2])
      .set{ norm_tumor_abra_bams_w_is_metrics }
    picard_mark_duplicates_w_mate_cigar(
      norm_tumor_abra_bams_w_is_metrics,
      params.somatic$samps_out_dir,
      params.shared_dir)
    index_feature_file(
      params.somatic$known_sites,
      params.somatic$indices_out_dir,
      params.shared_dir)
    base_recalibrator(
      picard_mark_duplicates_w_mate_cigar.out.alns,
      collective_idx_files,
      index_feature_file.out.ff_w_index,
      params.somatic$raw_alns_to_cadabra_somatic_vcfs$base_recalibrator_params,
      params.somatic$samps_out_dir,
      params.shared_dir)
    apply_bqsr(
      base_recalibrator.out.grps,
      collective_idx_files,
      params.somatic$raw_alns_to_cadabra_somatic_vcfs$apply_bqsr_params,
      params.somatic$samps_out_dir,
      params.shared_dir)
    samtools_index_bqsr(
      apply_bqsr.out.alns,
      params.somatic$samps_out_dir,
      params.shared_dir)
    // This should probably be its own workflow...
    manifest.filter{ it[5] =~ /[Nn]orm/ }.set{ norms }
    manifest.filter{ it[5] =~ /[Tt]umor/ }.set{ tumors }
    apply_bqsr.out.alns
      .join(samtools_index_bqsr.out.bais, by: [0, 1, 2])
      .set{ procd_bams_bais }
    procd_bams_bais.join(norms, by: [0, 1, 2])
      .map{ [it[0], it[1], it[2], it[3], it[4]] }
      .set{ norm_procd_bams_bais }
    procd_bams_bais.join(tumors, by: [0, 1, 2])
      .map{ [it[0], it[1], it[2], it[3], it[4]] }
      .set{ tumor_procd_bams_bais }
    norm_procd_bams_bais
      .join(tumor_procd_bams_bais, by: [0, 2])
      .set{ norm_tumor_procd_bams_bais }
    procd_alns_to_cadabra_somatic_vcfs(norm_tumor_procd_bams_bais, dna_ref)
  emit:
    vcfs = procd_alns_to_cadabra_somatic_vcfs.out.vcfs
}


workflow procd_alns_to_mutect2_somatic_vcfs {
// require
//   BAMS_AND_BAIS
//   REF_W_INDICES
  take:
    bams_and_bais
    ref_w_indices
    targets_bed
    mutect2_pon_vcf
    mutect2_af_vcf
    mutect2_parameters
    mutect2_suffix
    species
  main:
    htslib_tabix_pon(
      mutect2_pon_vcf)
    htslib_tabix_af(
      mutect2_af_vcf)
    gatk_mutect2_matched(
      bams_and_bais,
      ref_w_indices,
      targets_bed,
      mutect2_parameters,
      mutect2_suffix,
      htslib_tabix_pon.out.tabix_file,
      htslib_tabix_af.out.tabix_file,
      species)
  emit:
    vcfs = gatk_mutect2_matched.out.vcfs
    vcfs_w_stats = gatk_mutect2_matched.out.vcfs_w_stats
    f1r2_tar_gzs = gatk_mutect2_matched.out.f1r2_tar_gzs
}


workflow procd_alns_to_cadabra_somatic_vcfs {
// require:
//   BAMS_AND_BAIS
//   REF
  take:
    bams_and_bais
    ref
    abra2_cadabra_parameters
    abra2_cadabra_suffix
  main:
    abra2_cadabra(
      bams_and_bais,
      ref,
      abra2_cadabra_suffix,
      abra2_cadabra_parameters)
  emit:
    vcfs = abra2_cadabra.out.vcfs
}


workflow procd_alns_to_strelka2_somatic_vcfs {
  take:
    bams_bais
    ref
    ref_w_indices
    target_bed_w_tbi
    strelka2_parameters
    strelka2_suffix
  main:
    strelka2_somatic(
      bams_bais,
      ref_w_indices,
      target_bed_w_tbi,
      strelka2_parameters,
      strelka2_suffix)
  emit:
    snv_vcfs = strelka2_somatic.out.snv_vcfs
    indel_vcfs = strelka2_somatic.out.indel_vcfs
}


workflow bams_to_base_qual_recal_w_indices {
  take:
    alns
    known_sites
    collective_idx_files
    samtools_index_parameters
    gatk_index_feature_file_parameters
    targets_bed
    gatk_base_recalibrator_parameters
    gatk_apply_bqsr_parameters
  main:
    samtools_index_recal(
      alns,
      samtools_index_parameters)
    gatk_index_feature_file(
      known_sites,
      gatk_index_feature_file_parameters)
    gatk_base_recalibrator(
      samtools_index_recal.out.bams_and_bais,
      collective_idx_files,
      gatk_index_feature_file.out.ff_w_index,
      targets_bed,
      gatk_base_recalibrator_parameters)
    gatk_apply_bqsr(
      gatk_base_recalibrator.out.grps,
      collective_idx_files,
      targets_bed,
      gatk_apply_bqsr_parameters)
    samtools_index_bqsr(
      gatk_apply_bqsr.out.bams,
      samtools_index_parameters)
  emit:
    bams_and_bais = samtools_index_bqsr.out.bams_and_bais
    bams = gatk_apply_bqsr.out.bams
    bais = samtools_index_bqsr.out.bais
}


//Maybe this shouldn't be in somatic, but rather in a generic alignments tools module?
//I think this is identical is get_ref_idx_files
workflow make_ancillary_index_files {
//This creates index files useful for variant calling. Currently, this includes
//the samtools FASTA index file and Picard SequenceDictionary file.
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
