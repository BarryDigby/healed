include { whatshap } from '../whatshap/whatshap.nf'
include { whatshap_tumor } from '../whatshap/whatshap.nf'
include { whatshap_stats } from '../whatshap/whatshap.nf'
include { bcftools_view } from '../bcftools/bcftools.nf'
include { bcftools_reheader as bcftools_reheader_somatic } from '../bcftools/bcftools.nf'
include { bcftools_reheader as bcftools_reheader_germline } from '../bcftools/bcftools.nf'
include { bcftools_concat } from '../bcftools/bcftools.nf'
include { bcftools_index_somatic } from '../bcftools/bcftools.nf'
include { bcftools_index_somatic as bcftools_index_germline } from '../bcftools/bcftools.nf'
include { bcftools_index_somatic as bcftools_index_phased_germline } from '../bcftools/bcftools.nf'
include { bcftools_index_somatic as bcftools_index_tumor } from '../bcftools/bcftools.nf'
include { bcftools_index_somatic as bcftools_index_phased_tumor } from '../bcftools/bcftools.nf'
include { htslib_bgzip_somatic } from '../htslib/htslib.nf'
include { htslib_bgzip_somatic as htslib_bgzip_germline} from '../htslib/htslib.nf'
include { htslib_bgzip_somatic as htslib_bgzip_phased_germline} from '../htslib/htslib.nf'
include { htslib_bgzip_somatic as htslib_bgzip_tumor } from '../htslib/htslib.nf'
include { htslib_bgzip_somatic as htslib_bgzip_phased_tumor} from '../htslib/htslib.nf'

include { lenstools_get_expressed_transcripts_bed } from '../lenstools/lenstools.nf'
include { lenstools_get_fusion_transcripts_bed } from '../lenstools/lenstools.nf'
include { samtools_faidx_fetch_somatic_folder } from '../samtools/samtools.nf'
include { samtools_faidx_fetch } from '../samtools/samtools.nf'

include { lenstools_make_variant_specific_vcfs } from '../lenstools/lenstools.nf'
include { lenstools_make_variant_specific_tx_seqs } from '../lenstools/lenstools.nf'
include { starfusion_fusions_to_fusion_txs } from '../fusion/fusion.nf'

include { lenstools_make_fusion_specific_tx_seqs } from '../lenstools/lenstools.nf'

workflow tx_list_to_norm_and_tumor_cds_fastas {
  take:
    somatic_vcfs
    tx_lists
    gtf
    dna_ref
    normal_phased_vcfs
    tumor_phased_vcfs
    bcftools_index_phased_germline_parameters
    bcftools_index_phased_tumor_parameters
    lenstools_get_expressed_transcripts_bed_parameters
    samtools_faidx_fetch_somatic_folder_parameters
  main:
    htslib_bgzip_phased_germline(
      normal_phased_vcfs)
    bcftools_index_phased_germline(
      htslib_bgzip_phased_germline.out.bgzip_files,
      bcftools_index_phased_germline_parameters)
    htslib_bgzip_phased_tumor(
      tumor_phased_vcfs)
    bcftools_index_phased_tumor(
      htslib_bgzip_phased_tumor.out.bgzip_files,
      bcftools_index_phased_tumor_parameters)
    lenstools_get_expressed_transcripts_bed(
      tx_lists,
      gtf,
      lenstools_get_expressed_transcripts_bed_parameters)
    somatic_vcfs
      .join(normal_phased_vcfs, by: [0, 1, 2, 3])
      .join(tumor_phased_vcfs, by: [0, 1, 2, 3])
      .set{ somatic_w_norm_w_tumor }
    lenstools_make_variant_specific_vcfs(
      somatic_w_norm_w_tumor)
    samtools_faidx_fetch_somatic_folder(
      lenstools_get_expressed_transcripts_bed.out.expressed_transcripts_beds,
      dna_ref,
      'expressed_txs',
      samtools_faidx_fetch_somatic_folder_parameters)
    lenstools_make_variant_specific_vcfs.out.somatic_vcf_and_vcf_spec_vcfs
      .join(samtools_faidx_fetch_somatic_folder.out.fetched_fastas, by: [0, 1, 2, 3])
      .set{ var_vcfs_and_tx_seqs }
    lenstools_make_variant_specific_tx_seqs(
      var_vcfs_and_tx_seqs)
  emit:
    somatic_vcfs_w_var_tx_seqs = lenstools_make_variant_specific_tx_seqs.out.somatic_vcfs_w_var_tx_seqs
}

workflow create_phased_germline_tx_vcfs {
// bams
// vcf
// gtf
// tx list
  take:
    vcfs_and_csis
    bams_and_bais
    tx_lists
    ref_and_faidx
    gtf
    whatshap_parameters
    whatshap_stats_parameters
  main:
    vcfs_and_csis
      .join(bams_and_bais, by: [0, 1, 2, 3])
      .set{ vcfs_and_bams }
    whatshap(
      vcfs_and_bams,
      ref_and_faidx,
      whatshap_parameters)
    whatshap_stats(
      whatshap.out.phased_vcfs,
      whatshap_stats_parameters)
  emit:
    phased_vcfs = whatshap.out.phased_vcfs
    phased_stats = whatshap_stats.out.phased_stats
}


workflow create_phased_tumor_tx_vcfs {
// bams
// vcf
// gtf
// tx list
  take:
    germline_vcfs_and_csis
    somatic_vcfs_and_csis
    dna_bams_and_bais
    rna_bams_and_bais
    tx_lists
    ref_and_faidx
    gtf
    bcftools_view_parameters
    bcftools_reheader_somatic_parameters
    bcftools_index_somatic_parameters
    bcftools_reheader_germline_parameters
    bcftools_index_germline_parameters
    bcftools_concat_parameters
    bcftools_index_tumor_parameters
    whatshap_tumor_parameters
    whatshap_stats_parameters
  main:
    bcftools_view(
      somatic_vcfs_and_csis,
      'somatic_only',
      bcftools_view_parameters)
    htslib_bgzip_somatic(
      bcftools_view.out.viewed_vcfs)
    bcftools_reheader_somatic(
      htslib_bgzip_somatic.out.bgzip_files,
      bcftools_reheader_somatic_parameters)
    bcftools_index_somatic(
      bcftools_reheader_somatic.out.reheadered_vcfs,
      bcftools_index_somatic_parameters)
    bcftools_reheader_germline(
      germline_vcfs_and_csis.map{ [it[0], it[1], it[2], it[3], it[4]] },
      bcftools_reheader_germline_parameters)
    bcftools_index_germline(
      bcftools_reheader_germline.out.reheadered_vcfs,
      bcftools_index_germline_parameters)
    bcftools_index_somatic.out.vcfs_w_csis
      .join(bcftools_index_germline.out.vcfs_w_csis, by: [0, 1, 2, 3])
      .map{ [it[0], it[1], it[2], it[3], [it[4], it[6]], [it[5], it[7]]] }
      .set{ somatic_and_germline_vcfs_w_csis }
    bcftools_concat(
      somatic_and_germline_vcfs_w_csis,
      bcftools_concat_parameters, //'-a',
      'tumor')
    htslib_bgzip_tumor(
      bcftools_concat.out.concatd_vcfs)
    bcftools_index_tumor(
      htslib_bgzip_tumor.out.bgzip_files,
      bcftools_index_tumor_parameters)
    bcftools_index_tumor.out.vcfs_w_csis
      .join(dna_bams_and_bais, by: [0, 1, 2, 3])
      .set{ concatd_vcfs_and_dna_bams }
    concatd_vcfs_and_dna_bams.join(rna_bams_and_bais, by: 0)
      .map{ [it[0], it[1], it[2], it[8], it[3], it[4], it[5], it[6], it[7], it[10], it[11]] }
      .set{ concatd_vcfs_and_dna_bams_and_rna_bams }
    whatshap_tumor(
      concatd_vcfs_and_dna_bams_and_rna_bams,
      ref_and_faidx,
      whatshap_tumor_parameters) //'--ignore-read-groups'
    whatshap_stats(
      whatshap_tumor.out.phased_vcfs,
      whatshap_stats_parameters)
  emit:
    phased_vcfs = whatshap_tumor.out.phased_vcfs
    phased_stats = whatshap_stats.out.phased_stats
}

//process extract_transcripts_from_fusions {
//
//  input:
//  tuple val(pat_name), val(run), val(dataset), path(fusions)
//
//  output:
//  tuple val(path(name), val(run), val(dataset), path("*fusion_transcripts.txt"), emit: fusion_transcripts
//
//  script:
//  """
//  echo "Hello"
//  """
//}


workflow fusions_to_fusion_cds_fastas {
  take:
    fusions
    gtf
    dna_ref
    normal_phased_vcfs
    bcftools_index_phased_germline_parameters
    lenstools_get_fusion_transcripts_bed_parameters
    samtools_faidx_fetch_parameters
  main:
    htslib_bgzip_phased_germline(
      normal_phased_vcfs)
    bcftools_index_phased_germline(
      htslib_bgzip_phased_germline.out.bgzip_files,
      bcftools_index_phased_germline_parameters)
    bcftools_index_phased_germline.out.vcfs_w_csis
      .map{ [it[0], it[1], it[3], it[4], it[5]] }
      .set{ normal_phased_vcfs_trunc }
    fusions.combine(normal_phased_vcfs_trunc, by: [0, 2])
      .map{ [it[0], it[2], it[4], it[1], it[3], it[5], it[6]] }
      .set{ fusions_w_normal_phased_vcfs }
    starfusion_fusions_to_fusion_txs(
      fusions)
    lenstools_get_fusion_transcripts_bed(
      starfusion_fusions_to_fusion_txs.out.fusion_transcripts,
      gtf,
      lenstools_get_fusion_transcripts_bed_parameters)
    samtools_faidx_fetch(
      lenstools_get_fusion_transcripts_bed.out.fusion_transcripts_beds,
      dna_ref,
      'expressed_txs',
      samtools_faidx_fetch_parameters)
    samtools_faidx_fetch.out.fetched_fastas
      .combine(normal_phased_vcfs_trunc, by: [0, 2])
      .map{ [it[0], it[2], it[4], it[1], it[3], it[5], it[6]] }
      .set{ exon_fas_w_normal_phased_vcfs }
    lenstools_make_fusion_specific_tx_seqs(
      exon_fas_w_normal_phased_vcfs)
  emit:
      fusion_tx_exons = lenstools_make_fusion_specific_tx_seqs.out.var_tx_seqs
}
