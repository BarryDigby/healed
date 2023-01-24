#!/usr/bin/env nextflow

// Preprocessing
include { manifest_to_rna_procd_fqs } from '../preproc/preproc.nf'
include { manifest_to_dna_procd_fqs } from '../preproc/preproc.nf'
include { procd_fqs_to_procd_alns } from '../somatic/somatic.nf'

// Pre-alignment filtering
//include { filter_channel_by_manifest as filter_channel_by_manifest_norm_dnas } from '../utilities/utilities.nf'

// QC
include { initial_qc } from '../qc/qc.nf'
include { bams_to_rna_seq_metrics } from '../qc/qc.nf'
include { bams_to_rna_seq_metrics as realigned_bams_to_rna_seq_metrics } from '../qc/qc.nf'
include { vcfs_to_vcf_metrics } from '../qc/qc.nf'
include { samtools_stats } from  '../samtools/samtools.nf'
include { samtools_index } from  '../samtools/samtools.nf'
include { samtools_index as samtools_index_alt} from  '../samtools/samtools.nf'
include { samtools_sort as samtools_sort_alt} from  '../samtools/samtools.nf'
include { samtools_index as samtools_index_viral} from  '../samtools/samtools.nf'
include { samtools_index as samtools_index_prerealign } from  '../samtools/samtools.nf'
include { samtools_sort as samtools_sort_viral} from  '../samtools/samtools.nf'
include { samtools_coverage as samtools_coverage_viral} from  '../samtools/samtools.nf'
include { samtools_coverage as samtools_coverage_erv} from  '../samtools/samtools.nf'
include { bcftools_stats_somatic} from  '../bcftools/bcftools.nf'
include { bcftools_stats as bcftools_stats_germline } from  '../bcftools/bcftools.nf'
include { multiqc } from '../multiqc/multiqc.nf'

// Transcript quantification
include { procd_fqs_to_star_alns_salmon_counts } from '../rna_quant/rna_quant.nf'

// Transcript alignments (includes indel realignemnt)
include { procd_fqs_to_indel_realigned_star_alns } from '../rna_quant/rna_quant.nf'
include { abra2_rna } from '../abra2/abra2.nf'

// Viral expression
include { alns_to_virdetect } from '../viral/viral.nf'
include { unaligned_fqs_to_virdetect_cds_counts } from '../viral/viral.nf'
include { lenstools_filter_ervs_by_rna_coverage } from '../lenstools/lenstools.nf'
include { lenstools_filter_viruses_by_rna_coverage } from '../lenstools/lenstools.nf'

// Neosplice (splice variants, retained introns, etc.)
include { bams_to_neosplice } from '../rna_antigens/rna_antigens.nf'
include { neosplice_filter_and_summarize_variants } from '../neosplice/neosplice.nf'
include { lenstools_get_splice_peptide_read_count } from '../lenstools/lenstools.nf'
include { antigen_garnish_foreignness_rna } from '../antigen.garnish/antigen.garnish.nf'
include { antigen_garnish_dissimilarity_rna } from '../antigen.garnish/antigen.garnish.nf'
include { netmhcstabpan_rna } from '../netmhcstabpan/netmhcstabpan.nf'

// MHC typing
include { procd_fqs_to_hlaprofiler_calls } from '../immuno/immuno.nf'
include { procd_fqs_to_optitype_calls } from '../immuno/immuno.nf'
include { hlaprofiler_alleles_to_netmhcpan_alleles } from '../immuno/immuno.nf'
include { optitype_alleles_to_netmhcpan_alleles } from '../immuno/immuno.nf'
include { extract_alleles_from_manifest } from '../immuno/immuno.nf'
include { user_provided_alleles_to_netmhcpan_alleles } from '../immuno/immuno.nf'

// Fusion prediction
//include { junctions_to_starfusion_fusions } from '../fusion/fusion.nf'
include { procd_fqs_to_starfusion_fusions } from '../fusion/fusion.nf'

// Somatic variant calling
include { procd_alns_to_mutect2_somatic_vcfs } from '../somatic/somatic.nf'
include { procd_alns_to_strelka2_somatic_vcfs } from '../somatic/somatic.nf'
include { procd_alns_to_cadabra_somatic_vcfs } from '../somatic/somatic.nf'

// Germline variant calling
include { procd_alns_to_deepvariant_germline_vcfs } from '../germline/germline.nf'

// VCF normalization and intersection
include { bcftools_norm as bcftools_norm_cadabra } from '../bcftools/bcftools.nf'
include { bcftools_norm as bcftools_norm_mutect2 } from '../bcftools/bcftools.nf'
include { bcftools_norm as bcftools_norm_strelka2 } from '../bcftools/bcftools.nf'
include { bcftools_isec } from '../bcftools/bcftools.nf'
include { bedtools_intersect_vcfs } from '../bedtools/bedtools.nf'
include { bcftools_index } from '../bcftools/bcftools.nf'
include { bcftools_index_somatic as bcftools_index_isec} from '../bcftools/bcftools.nf'
include { bcftools_index_somatic } from '../bcftools/bcftools.nf'

// Variant annotation
include { snpeff_ann } from '../snpeff/snpeff.nf'
include { snpeff_ann_germline } from '../snpeff/snpeff.nf'

// Variant filtering
include { snpsift_filter as snpsift_filter_snvs } from '../snpeff/snpeff.nf'
include { snpsift_filter as snpsift_filter_indels} from '../snpeff/snpeff.nf'
include { bcftools_filter as bcftools_filter_mutect2 } from '../bcftools/bcftools.nf'
include { bcftools_filter as bcftools_filter_cadabra } from '../bcftools/bcftools.nf'
include { bcftools_filter as bcftools_filter_snvs_strelka2 } from '../bcftools/bcftools.nf'
include { bcftools_filter as bcftools_filter_indels_strelka2 } from '../bcftools/bcftools.nf'
include { bcftools_filter_germline } from '../bcftools/bcftools.nf'
include { gatk_filter_mutect_calls } from '../gatk4/gatk4.nf'
include { gatk_get_pileup_summaries } from '../gatk4/gatk4.nf'
include { gatk_calculate_contamination } from '../gatk4/gatk4.nf'
include { gatk_index_feature_file } from '../gatk4/gatk4.nf'
include { gatk_learn_read_orientation_model } from '../gatk4/gatk4.nf'

// Variant misc.
include { htslib_bgzip } from '../htslib/htslib.nf'
include { htslib_bgzip_ref } from '../htslib/htslib.nf'
include { htslib_tabix } from '../htslib/htslib.nf'
include { htslib_bgzip_somatic } from '../htslib/htslib.nf'
include { htslib_bgzip_somatic as htslib_bgzip_somatic_isec } from '../htslib/htslib.nf'
include { htslib_bgzip_somatic as htslib_bgzip_somatic_isec_annot } from '../htslib/htslib.nf'
include { filter_channel_by_manifest } from '../utilities/utilities.nf'
include { bcftools_merge } from '../bcftools/bcftools.nf'
include { create_phased_germline_tx_vcfs } from '../seq_variation/seq_variation.nf'
include { create_phased_tumor_tx_vcfs } from '../seq_variation/seq_variation.nf'

// Tumor ploidy estimation and CNV quantification
include { procd_alns_to_sequenza } from '../onco/onco.nf'

// Potential antigen variant filtering (e.g. proximity, expression)
include { lenstools_filter_expressed_variants } from '../lenstools/lenstools.nf'
include { lenstools_filter_isolated_variants } from '../lenstools/lenstools.nf'
include { lenstools_filter_expressed_ervs } from '../lenstools/lenstools.nf'
include { lenstools_filter_expressed_ervs_without_control } from '../lenstools/lenstools.nf'
include { lenstools_filter_expressed_viruses } from '../lenstools/lenstools.nf'
include { lenstools_get_viral_cds_expression } from '../lenstools/lenstools.nf'
include { lenstools_filter_expressed_self_genes } from '../lenstools/lenstools.nf'
include { lenstools_make_pyclonevi_inputs } from '../lenstools/lenstools.nf'
include { lenstools_add_splice_metadata } from '../lenstools/lenstools.nf'

// Tumor antigen prediction
include { filtered_snvs_to_class1_neos } from '../neos/neos.nf'
include { filtered_indels_to_class1_neos } from '../neos/neos.nf'
include { filtered_ervs_to_class1_neos } from '../neos/neos.nf'
include { filtered_viruses_to_class1_neos } from '../neos/neos.nf'
include { filtered_selfs_to_class1_neos } from '../neos/neos.nf'
include { filtered_fusions_to_class1_neos } from '../neos/neos.nf'

// Cancer cell fraction
include { pyclonevi_fit } from '../pyclone-vi/pyclone-vi.nf'
include { pyclonevi_write_results_file } from '../pyclone-vi/pyclone-vi.nf'
include { lenstools_split_pyclonevi_outputs } from '../lenstools/lenstools.nf'

// Utilities
include { combine_sample_files } from '../utilities/utilities.nf'
include { combine_patient_samples } from '../utilities/utilities.nf'
include { combine_rel_patient_samples } from '../utilities/utilities.nf'
include { combine_rel_patient_samples as combine_rel_patient_samples_pileups } from '../utilities/utilities.nf'

// Misc.
include { samtools_faidx } from '../samtools/samtools.nf'
include { tx_list_to_norm_and_tumor_cds_fastas as tx_list_to_snv_norm_and_tumor_cds_fastas } from '../seq_variation/seq_variation.nf'
include { tx_list_to_norm_and_tumor_cds_fastas as tx_list_to_indel_norm_and_tumor_cds_fastas } from '../seq_variation/seq_variation.nf'
include { lenstools_make_lens_report } from '../lenstools/lenstools.nf'
include { lenstools_make_lens_report as lenstools_make_lens_reports_rna_only } from '../lenstools/lenstools.nf'
include { lenstools_make_lens_bed } from '../lenstools/lenstools.nf'
include { lenstools_add_tcga_data as lenstools_add_tcga_data_snvs } from '../lenstools/lenstools.nf'
include { lenstools_add_tcga_data as lenstools_add_tcga_data_indels } from '../lenstools/lenstools.nf'
include { make_ancillary_index_files } from '../rna_quant/rna_quant.nf'

// Visualization
include { igv_snapshot_automator } from '../igv_snapshot_automator/igv_snapshot_automator.nf'
include { bam_subsetter } from '../samtools/samtools.nf'

workflow manifest_to_lens {
// require:
//   MANIFEST
  take:
    manifest

  main:
    // Downsampled QC
//    initial_qc(
//      manifest)

    // Preprocessing
    manifest_to_rna_procd_fqs(
      manifest,
      params.manifest_to_lens$manifest_to_rna_procd_fqs$trim_galore_parameters)
    manifest_to_dna_procd_fqs(
      manifest,
      params.manifest_to_lens$manifest_to_rna_procd_fqs$trim_galore_parameters)

    procd_fqs_to_procd_alns(
        manifest_to_dna_procd_fqs.out.procd_fqs,
        manifest,
        params.manifest_to_lens$procd_fqs_to_procd_alns$dna_ref,
        params.manifest_to_lens$procd_fqs_to_procd_alns$bwa_index_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$bwa_mem_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$targets_bed,
        params.manifest_to_lens$procd_fqs_to_procd_alns$picard_mark_duplicates_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$known_sites_ref,
        params.manifest_to_lens$procd_fqs_to_procd_alns$samtools_index_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$gatk_index_feature_file_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$gatk_base_recalibrator_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$gatk_apply_bqsr_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$samtools_faidx_parameters,
        params.manifest_to_lens$procd_fqs_to_procd_alns$abra2_parameters)


/* RNA-Sequencing processing */

    // Transcript quantification
    procd_fqs_to_star_alns_salmon_counts(
      params.manifest_to_lens$procd_fqs_to_star_alns_salmon_counts$star_salmon_ref,
      params.manifest_to_lens$procd_fqs_to_star_alns_salmon_counts$gtf,
      manifest_to_rna_procd_fqs.out.procd_fqs,
      params.manifest_to_lens$procd_fqs_to_star_alns$star_index_parameters,
      params.manifest_to_lens$procd_fqs_to_star_alns$star_map_parameters,
      params.manifest_to_lens$procd_fqs_to_star_alns$star_alt_capture,
      params.manifest_to_lens$procd_fqs_to_star_alns_salmon_counts$star_aln_quant_parameters)
    bams_to_rna_seq_metrics(
      params.manifest_to_lens$bams_to_rna_seq_metrics$star_salmon_ref,
      procd_fqs_to_star_alns_salmon_counts.out.alns,
      params.manifest_to_lens$bams_to_rna_seq_metrics$gtf,
      params.manifest_to_lens$bams_to_rna_seq_metrics$picard_collect_rna_seq_metrics_parameters,
      params.manifest_to_lens$bams_to_rna_seq_metrics$picard_strand_specificity,
      params.manifest_to_lens$bams_to_rna_seq_metrics$samtools_faidx_parameters)

/* MHC typing */

    // MHC typing
    if( params.mhc_alleles_manifest ) {
      extract_alleles_from_manifest(
        params.mhc_alleles_manifest)
      extract_alleles_from_manifest.out.patient_alleles.filter{ it[3] =~ /NA|Null|null/ }.set{ pats_missing_alleles }
      extract_alleles_from_manifest.out.patient_alleles.filter{ !(it[3] =~ /NA|Null|null/) }.set{ pats_with_alleles }
      if ( params.hla_call_method =~ /RNA|RNA-Seq|rna|rna-seq/ ) {
        procd_fqs_to_hlaprofiler_calls(
          pats_missing_alleles.filter{ it[1] =~ 'ar-' },
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_read_count,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_seed,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_sample_suffix,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_parameters,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$trim_galore_hlap_parameters,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$hlaprofiler_parameter)
        hlaprofiler_alleles_to_netmhcpan_alleles(
          procd_fqs_to_hlaprofiler_calls.out.calls)
          user_provided_alleles_to_netmhcpan_alleles(
          pats_with_alleles)
        hlaprofiler_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles
          .concat(user_provided_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles)
          .set{ netmhcpan_alleles }
      } else if ( params.hla_call_method =~ /WES|wes|WXS|wxs|DNA|dna/ ) {
        procd_fqs_to_optitype_calls(
          pats_missing_alleles.filter{ it[1] =~ 'nd-' })
        optitype_alleles_to_netmhcpan_alleles(
          procd_fqs_to_optitype_calls.out.calls)
          user_provided_alleles_to_netmhcpan_alleles(
          pats_with_alleles)
        optitype_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles
          .concat(user_provided_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles)
          .set{ netmhcpan_alleles }
      } // This comment is to keep this '}' from upsetting RAFT.
    } else {
      if ( params.hla_call_method =~ /RNA|RNA-Seq|rna|rna-seq/ ) {
        procd_fqs_to_hlaprofiler_calls(
          manifest_to_rna_procd_fqs.out.procd_fqs.filter{ it[1] =~ 'ar-' },
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_read_count,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_seed,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_sample_suffix,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$seqtk_parameters,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$trim_galore_hlap_parameters,
          params.manifest_to_lens$procd_fqs_to_hlaprofiler_calls$hlaprofiler_parameter)
        hlaprofiler_alleles_to_netmhcpan_alleles(
          procd_fqs_to_hlaprofiler_calls.out.calls)
        hlaprofiler_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles
          .set{ netmhcpan_alleles }
      } else if ( params.hla_call_method =~ /WES|wes|WXS|wxs|DNA|dna/ ) {
        procd_fqs_to_optitype_calls(
          manifest_to_dna_procd_fqs.out.procd_fqs.filter{ it[1] =~ 'nd-' })
        optitype_alleles_to_netmhcpan_alleles(
          procd_fqs_to_optitype_calls.out.calls)
        optitype_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles
          .set{ netmhcpan_alleles }
      } // This comment is to keep this '}' from upsetting RAFT.
    } // This comment is to keep this '}' from upsetting RAFT.

/* ERV/Viral expression */

    // Tumor Viral expression
    alns_to_virdetect(
      procd_fqs_to_star_alns_salmon_counts.out.alns.filter{ it[1] =~ 'ar-' },
      params.manifest_to_lens$alns_to_virdetect$virdetect_viral_ref)
    unaligned_fqs_to_virdetect_cds_counts(
      alns_to_virdetect.out.unaligned_fqs,
      params.manifest_to_lens$unaligned_fqs_to_virdetect_cds_counts$viral_cds_ref)


/* Fusion prediction */

    // Fusion detection
    procd_fqs_to_starfusion_fusions(
      params.manifest_to_lens$ctat_ref,
      params.manifest_to_lens$starfusion_parameters,
      manifest_to_rna_procd_fqs.out.procd_fqs.filter({it[1] =~ 'ar-' }))



/* Variant calling */
    // Indexing FASTA (for bcftools_norm)
    samtools_faidx(
      params.manifest_to_lens$samtools_faidx$dna_ref,
      params.manifest_to_lens$samtools_faidx$samtools_faidx_parameters)

    // Somatic variant calling
    combine_sample_files(
      procd_fqs_to_procd_alns.out.filt_bams,
      procd_fqs_to_procd_alns.out.filt_bais)
    combine_rel_patient_samples(
      combine_sample_files.out.combined_set)

    // Tumor ploidy, purity, and copy number
    procd_alns_to_sequenza(
      combine_rel_patient_samples.out.combined_rel_set,
      params.manifest_to_lens$procd_alns_to_sequenza$dna_ref,
      params.manifest_to_lens$procd_alns_to_sequenza$targets_bed,
      params.manifest_to_lens$procd_alns_to_sequenza$sequenza_gc_wiggle_window_size,
      params.manifest_to_lens$procd_alns_to_sequenza$sequenza_gc_wiggle_parameters,
      params.manifest_to_lens$procd_alns_to_sequenza$sequenza_seqz_binning_window_size)

    // MuTect2
    procd_alns_to_mutect2_somatic_vcfs(
      combine_rel_patient_samples.out.combined_rel_set,
      procd_fqs_to_procd_alns.out.anc_idx_files,
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$targets_bed,
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$mutect2_pon_vcf,
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$mutect2_af_vcf,
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$mutect2_parameters,
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$mutect2_suffix,
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$species)
    if( params.species =~ /[Hh]uman|hs|HS|[Hh]omo/ ) {
        chrom_count = 25
    } else if( params.species =~ /[Mm]ouse|mm|MM|[Mm]us/ ) {
        chrom_count = 21
    }  // Keep RAFT happy
    gatk_learn_read_orientation_model(
      procd_alns_to_mutect2_somatic_vcfs.out.f1r2_tar_gzs.groupTuple(by: [0,1,2,3], size:chrom_count),
      params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$gatk_learn_read_orientation_model_parameters)

    if( params.sites_file != "") {
      gatk_index_feature_file(
        params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$sites_file,
        params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$gatk_index_feature_file_parameters)
      gatk_get_pileup_summaries(
        combine_sample_files.out.combined_set,
        gatk_index_feature_file.out.ff_w_index,
        params.manifest_to_lens$procd_alns_to_mutect2_somatic_vcfs$get_pileup_summaries_parameters)
      combine_rel_patient_samples_pileups(
        gatk_get_pileup_summaries.out.pileups_tables)
      gatk_calculate_contamination(
        combine_rel_patient_samples_pileups.out.combined_rel_set.map{ [it[0], it[1], it[2], it[3], it[5], it[6]] },
        '')
      procd_alns_to_mutect2_somatic_vcfs.out.vcfs_w_stats
        .join(gatk_calculate_contamination.out.contamination_tables, by: [0, 1, 2, 3])
        .join(gatk_learn_read_orientation_model.out.artifact_priors, by: [0, 1, 2, 3])
        .set{ vcfs_w_stats_w_contam_w_arti_priors }
    } else {
      procd_alns_to_mutect2_somatic_vcfs.out.vcfs_w_stats
        .join(gatk_learn_read_orientation_model.out.artifact_priors, by: [0, 1, 2, 3])
        .map{ [it[0], it[1], it[2], it[3], it[4], it[5], params.dummy_file, it[6]] }
        .set{ vcfs_w_stats_w_contam_w_arti_priors }
    } // Keep RAFT happy
    gatk_filter_mutect_calls(
      vcfs_w_stats_w_contam_w_arti_priors,
      procd_fqs_to_procd_alns.out.anc_idx_files,
      params.manifest_to_lens$gatk_filter_mutect_calls$gatk_filter_mutect_calls_parameters,
      params.manifest_to_lens$gatk_filter_mutect_calls$gatk_filter_mutect_calls_suffix)
    bcftools_filter_mutect2(
      gatk_filter_mutect_calls.out.filtd_vcfs,
      params.manifest_to_lens$bcftools_filter_mutect2$bcftools_filter_parameters)
    bcftools_norm_mutect2(
      bcftools_filter_mutect2.out.filtd_vcfs,
      samtools_faidx.out.faidx_file,
      params.manifest_to_lens$bcftools_norm_mutect2$bcftools_norm_parameters)

    // ABRA2 CADABRA
    procd_alns_to_cadabra_somatic_vcfs(
      combine_rel_patient_samples.out.combined_rel_set,
      params.manifest_to_lens$procd_alns_to_cadabra_somatic_vcfs$dna_ref,
      params.manifest_to_lens$procd_alns_to_cadabra_somatic_vcfs$abra2_cadabra_parameters,
      params.manifest_to_lens$procd_alns_to_cadabra_somatic_vcfs$abra2_cadabra_suffix)
    bcftools_filter_cadabra(
      procd_alns_to_cadabra_somatic_vcfs.out.vcfs,
      params.manifest_to_lens$bcftools_filter_cadabra$bcftools_filter_parameters)
    bcftools_norm_cadabra(
      bcftools_filter_cadabra.out.filtd_vcfs,
      samtools_faidx.out.faidx_file,
      params.manifest_to_lens$bcftools_norm_cadabra$bcftools_norm_parameters)

    // Strelka2
    htslib_bgzip_ref(
      params.manifest_to_lens$procd_alns_to_strelka2_somatic_vcfs$targets_bed)
    htslib_tabix(
      htslib_bgzip_ref.out.bgzip_files)
    procd_alns_to_strelka2_somatic_vcfs(
      combine_rel_patient_samples.out.combined_rel_set,
      params.manifest_to_lens$procd_alns_to_strelka2_somatic_vcfs$dna_ref,
      procd_fqs_to_procd_alns.out.anc_idx_files,
      htslib_tabix.out.tabix_file,
      params.manifest_to_lens$procd_alns_to_strelka2_somatic_vcfs$strelka2_parameters,
      params.manifest_to_lens$procd_alns_to_strelka2_somatic_vcfs$strelka2_suffix)
    bcftools_filter_snvs_strelka2(
      procd_alns_to_strelka2_somatic_vcfs.out.snv_vcfs,
      params.manifest_to_lens$bcftools_filter_snvs_strelka2$bcftools_filter_parameters)
    bcftools_filter_indels_strelka2(
      procd_alns_to_strelka2_somatic_vcfs.out.indel_vcfs,
      params.manifest_to_lens$bcftools$bcftools_filter_strelka2$bcftools_filter_parameters)
    bcftools_norm_strelka2(
      bcftools_filter_indels_strelka2.out.filtd_vcfs,
      samtools_faidx.out.faidx_file,
      params.manifest_to_lens$bcftools_norm_strelka2$bcftools_norm_parameters)

/* Somatic variants intersection */
    // Combining raw somatic VCFs for running stats
    gatk_filter_mutect_calls.out.filtd_vcfs
      .concat(procd_alns_to_cadabra_somatic_vcfs.out.vcfs,
         procd_alns_to_strelka2_somatic_vcfs.out.snv_vcfs,
         procd_alns_to_strelka2_somatic_vcfs.out.indel_vcfs)
      .set{ all_raw_somatic_vcfs }

    // Combining for indexing prior to isec
    bcftools_norm_mutect2.out.normd_vcfs
      .mix(bcftools_norm_cadabra.out.normd_vcfs,
         bcftools_filter_snvs_strelka2.out.filtd_vcfs,
         bcftools_norm_strelka2.out.normd_vcfs)
      .set{ all_somatic_vcfs }


    // Compressing variants
    htslib_bgzip_somatic(
      all_somatic_vcfs)

    // Indexing variants
    bcftools_index_somatic(
      htslib_bgzip_somatic.out.bgzip_files,
      '')

    // Collecting VCF metrics
//    vcfs_to_vcf_metrics(
//      bcftools_index_somatic.out.vcfs_w_csis,
//      params.manifest_to_lens$vcfs_to_vcf_metrics$dbsnp_ref,
//      params.manifest_to_lens$vcfs_to_vcf_metrics$picard_collect_vcf_metrics_parameters)

    // Group by sample
    bcftools_index_somatic.out.vcfs_w_csis.groupTuple(by: [0, 1, 2, 3], size: 4).set{ multi_vcfs }

    // Intersection for each sample
    bedtools_intersect_vcfs(
      multi_vcfs.map{ [it[0], it[1], it[2], it[3], it[4]] }) //Stripping off tabix indices

    htslib_bgzip_somatic_isec(
      bedtools_intersect_vcfs.out.intersect_vcfs)

    // Somatic variant annotation
    snpeff_ann(
      htslib_bgzip_somatic_isec.out.bgzip_files,
      params.manifest_to_lens$snpeff_ann$snpeff_ann_ref)

    htslib_bgzip_somatic_isec_annot(
      snpeff_ann.out.annot_vcfs)

    // Indexing somatic intersection VCF
    bcftools_index_isec(
      htslib_bgzip_somatic_isec_annot.out.bgzip_files,
      '-t')


    // Indexing RNA BAMs
    samtools_index_prerealign(
      procd_fqs_to_star_alns_salmon_counts.out.alns,
      params.manifest_to_lens$samtools_index$samtools_index_parameters)
    samtools_index_prerealign.out.bams_and_bais
      .join(procd_fqs_to_star_alns_salmon_counts.out.standard_junctions, by: [0, 1, 2])
      .map{ [it[0], it[2], it[1], it[3], it[4], it[5]] }
      .join(bedtools_intersect_vcfs.out.intersect_vcfs.map{ [it[0], it[3], it[4]] }, by: [0, 1])
      .map{ [it[0], it[2], it[1], it[3], it[4], it[5], it[6]] }
      .set{ bams_bais_junctions_vcfs }

    make_ancillary_index_files(
      params.star_salmon_ref)

    abra2_rna(
      bams_bais_junctions_vcfs,
      make_ancillary_index_files.out.collective_idx_files,
      params.gtf,
      params.targets_bed,
      '--dist 500000 --sua',
      'tmp_dir')


    samtools_index(
      abra2_rna.out.abra_bams
        .concat(procd_fqs_to_star_alns_salmon_counts.out.alns.filter{ it[1] =~ 'nr-' }),
      params.manifest_to_lens$samtools_index$samtools_index_parameters)


//* Splice variants detection */


    // Neosplice
    bams_to_neosplice(
      manifest.filter{ it[4] =~ /RNA/ },
//      samtools_index.out.bams_and_bais,
      samtools_index_prerealign.out.bams_and_bais,
      params.manifest_to_lens$bams_to_neosplice$neosplice_ref,
      params.manifest_to_lens$bams_to_neosplice$gff,
      netmhcpan_alleles,
      params.manifest_to_lens$bams_to_neosplice$neosplice_augmented_splice_graph_build_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_get_max_kmer_length_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_convert_bams_to_fasta_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_get_splice_junctions_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_msbwt_is_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_convert_bwt_format_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_kmer_search_bwt_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_search_bam_parameters,
      params.manifest_to_lens$bams_to_neosplice$samtools_sort_parameters,
      params.manifest_to_lens$bams_to_neosplice$samtools_index_parameters,
      params.manifest_to_lens$bams_to_neosplice$neosplice_kmer_graph_inference_parameters)

/* Calculating CFF from somatic variant calls */

    procd_alns_to_sequenza.out.segments_and_solutions.map{ [it[0], it[2], it[3], it[1], it[4], it[5]] }.set{ reorged_segments_and_solutions }

    //Replace with intersection VCFs
    bedtools_intersect_vcfs.out.intersect_vcfs
      .join(bcftools_filter_mutect2.out.filtd_vcfs, by: [0, 1, 2, 3])
      .join(procd_alns_to_sequenza.out.segments_and_solutions, by: [0, 1, 2, 3])
      .set{ joint_pcvi }


    lenstools_make_pyclonevi_inputs(
      joint_pcvi,
      params.manifest_to_lens$lenstools_make_pyclonevi_inputs_parameters)

    pyclonevi_fit(
      lenstools_make_pyclonevi_inputs.out.pcvi_inputs,
      params.manifest_to_lens$pyclonevi_fit$pyclonevi_fit_parameters)

    pyclonevi_write_results_file(
      pyclonevi_fit.out.pcvi_tmps,
      params.manifest_to_lens$pyclonevi_write_results_file$pyclonevi_write_results_file_parameters)

/* Germline variant calling */

    // Filtering for normal samples
    filter_channel_by_manifest(
      combine_sample_files.out.combined_set, //From variant calling set.
      '[Nn]orm',
      4,
      manifest)

    // Calling germline variants
    procd_alns_to_deepvariant_germline_vcfs(
      procd_fqs_to_procd_alns.out.filt_bams_and_bais.filter{ it[1].split('-rel-')[0] == it[1].split('-rel-')[1].split('-to-')[0] },
      params.manifest_to_lens$procd_alns_to_deepvariant_germline_vcfs$dna_ref,
      params.manifest_to_lens$procd_alns_to_deepvariant_germline_vcfs$targets_bed,
      params.manifest_to_lens$procd_alns_to_deepvariant_germline_vcfs$deepvariant_model_type,
      params.manifest_to_lens$procd_alns_to_deepvariant_germline_vcfs$deepvariant_parameters,
      params.manifest_to_lens$procd_alns_to_deepvariant_germline_vcfs$deepvariant_suffix,
      params.manifest_to_lens$procd_alns_to_deepvariant_germline_vcfs$samtools_faidx_parameters)

    // Germline variant filtering
    bcftools_filter_germline(
      procd_alns_to_deepvariant_germline_vcfs.out.germline_vcfs,
      params.manifest_to_lens$bcftools_filter_germline$bcftools_filter_parameters)


    // Germline variant compressing
    htslib_bgzip(
      bcftools_filter_germline.out.filtd_vcfs)

    // Germline variant annotation
    snpeff_ann_germline(
      htslib_bgzip.out.bgzip_files,
      params.manifest_to_lens$snpeff_ann$snpeff_ann_ref)

    // Germline variant indexing
    bcftools_index(
      htslib_bgzip.out.bgzip_files,
      params.manifest_to_lens$bcftools_index$bcftools_index_parameters)


/* Merge germline and somatic variants */
    // Expanding for joining...
    bcftools_index.out.vcfs_w_csis
      .map{ [it[0], it[1].split('-rel-')[1].split('-to-')[0], it[1].split('-rel-')[1].split('-to-')[1], it[2], it[3], it[4]] }
      .set{ bcftools_index_remapped }


/* Tumor antigen prediction */

    // Variant-based neoantigen prediction
    bcftools_index_isec.out.vcfs_w_csis
      .join(bcftools_index_remapped, by: [0, 1, 2, 3])
      .map{ [it[0], it[1], it[2], it[3], it[4], it[6]] }
      .set{ germline_and_somatic_joint_vcfs }

    // Filter for proximal variants
    lenstools_filter_isolated_variants(
      germline_and_somatic_joint_vcfs,
      params.manifest_to_lens$lenstools_filter_isolated_variants_parameters)

    // Filter for variant-associated transcript expression
    procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' }.set{ just_tumor_quants }
    lenstools_filter_isolated_variants.out.isolated_vcfs
      .combine(just_tumor_quants, by: [0])
      .map{ [it[0], it[1], it[2], it[3], it[4], it[7]] }
      .set{ annot_vars_and_exp }
    lenstools_filter_expressed_variants(
      annot_vars_and_exp,
      params.manifest_to_lens$lenstools_filter_expressed_variants_parameters)


    create_phased_tumor_tx_vcfs(
      bcftools_index_remapped,
      bcftools_index_isec.out.vcfs_w_csis,
      procd_fqs_to_procd_alns.out.filt_bams_and_bais.filter{ it[1].split('-rel-')[0] == it[1].split('-rel-')[1].split('-to-')[1] }.map{ [it[0], it[1].split('-rel-')[1].split('-to-')[0], it[1].split('-rel-')[1].split('-to-')[1], it[2], it[3], it[4]] },
      samtools_index.out.bams_and_bais,
      lenstools_filter_expressed_variants.out.somatic_transcripts,
      samtools_faidx.out.faidx_file,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$gtf,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_view_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_reheader_somatic_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_index_somatic_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_reheader_germline_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_index_germline_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_concat_parameters, // -a,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$bcftools_index_tumor_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$whatshap_tumor_parameters,
      params.manifest_to_lens$created_phased_tumor_tx_vcfs$whatshap_stats_parameters)

    create_phased_germline_tx_vcfs(
      bcftools_index_remapped,
      procd_fqs_to_procd_alns.out.filt_bams_and_bais
        .filter{ it[1].split('-rel-')[0] == it[1].split('-rel-')[1].split('-to-')[0] }
        .map{ [it[0], it[1].split('-rel-')[1].split('-to-')[0], it[1].split('-rel-')[1].split('-to-')[1], it[2], it[3], it[4]] },
      lenstools_filter_expressed_variants.out.somatic_transcripts,
      samtools_faidx.out.faidx_file,
      params.manifest_to_lens$created_phased_germline_tx_vcfs$gtf,
      params.manifest_to_lens$created_phased_germline_tx_vcfs$whatshap_parameters,
      params.manifest_to_lens$created_phased_germline_tx_vcfs$whatshap_stats_parameters)

    // Somatic SNV selection
    snpsift_filter_snvs(
      lenstools_filter_expressed_variants.out.expressed_vcfs,
      params.manifest_to_lens$snpsift_filter_snvs$snpsift_snv_filter_parameters,
      "sfilt.snvs")

    tx_list_to_snv_norm_and_tumor_cds_fastas(
      snpsift_filter_snvs.out.filt_vcfs,
      lenstools_filter_expressed_variants.out.somatic_transcripts,
      params.manifest_to_lens$tx_list_to_snv_norm_and_tumor_cds_fastas$gtf,
      params.manifest_to_lens$tx_list_to_snv_norm_and_tumor_cds_fastas$dna_ref,
      create_phased_germline_tx_vcfs.out.phased_vcfs,
      create_phased_tumor_tx_vcfs.out.phased_vcfs,
      params.manifest_to_lens$tx_list_to_snv_norm_and_tumor_cds_fastas$bcftools_index_phased_germline_parameters,
      params.manifest_to_lens$tx_list_to_snv_norm_and_tumor_cds_fastas$bcftools_index_phased_tumor_parameters,
      params.manifest_to_lens$tx_list_to_snv_norm_and_tumor_cds_fastas$lenstools_get_expressed_transcripts_bed_parameters,
      params.manifest_to_lens$tx_list_to_snv_norm_and_tumor_cds_fastas$samtools_faidx_fetch_somatic_folder_parameters)

    samtools_sort_alt(
      procd_fqs_to_star_alns_salmon_counts.out.alt_alns.filter{ it[1] =~ 'ar-' },
      '')

    samtools_index_alt(
      samtools_sort_alt.out.bams,
      '')

    // SNV antigen prediction
    filtered_snvs_to_class1_neos(
      snpsift_filter_snvs.out.filt_vcfs,
      tx_list_to_snv_norm_and_tumor_cds_fastas.out.somatic_vcfs_w_var_tx_seqs,
      lenstools_filter_expressed_variants.out.somatic_transcripts,
      netmhcpan_alleles,
      pyclonevi_write_results_file.out.pcvi_results,
      procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
      samtools_index_alt.out.bams_and_bais.filter{ it[1] =~ 'ar-' },
      params.manifest_to_lens$filtered_snvs_to_class1_neos$gtf,
      params.manifest_to_lens$filtered_snvs_to_class1_neos$netmhcpan_parameters,
      params.manifest_to_lens$filtered_snvs_to_class1_neos$netmhcstabpan_parameters,
      params.manifest_to_lens$filtered_snvs_to_class1_neos$lenstools_add_snv_metadata_parameters,
      params.manifest_to_lens$filtered_snvs_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_snvs_to_class1_neos$species,
      params.manifest_to_lens$filtered_snvs_to_class1_neos$protein_ref)


    // Somatic InDel selection
    snpsift_filter_indels(
      lenstools_filter_expressed_variants.out.expressed_vcfs,
      params.manifest_to_lens$snpeff$snpsift_indel_filter_parameters,
      "sfilt.indels")

    tx_list_to_indel_norm_and_tumor_cds_fastas(
      snpsift_filter_indels.out.filt_vcfs,
      lenstools_filter_expressed_variants.out.somatic_transcripts,
      params.manifest_to_lens$tx_list_to_indel_norm_and_tumor_cds_fastas$gtf,
      params.manifest_to_lens$tx_list_to_indel_norm_and_tumor_cds_fastas$dna_ref,
      create_phased_germline_tx_vcfs.out.phased_vcfs,
      create_phased_tumor_tx_vcfs.out.phased_vcfs,
      params.manifest_to_lens$tx_list_to_indel_norm_and_tumor_cds_fastas$bcftools_index_phased_germline_parameters,
      params.manifest_to_lens$tx_list_to_indel_norm_and_tumor_cds_fastas$bcftools_index_phased_tumor_parameters,
      params.manifest_to_lens$tx_list_to_indel_norm_and_tumor_cds_fastas$lenstools_get_expressed_transcripts_bed_parameters,
      params.manifest_to_lens$tx_list_to_indel_norm_and_tumor_cds_fastas$samtools_faidx_fetch_somatic_folder_parameters)

    // InDel antigen prediction
    filtered_indels_to_class1_neos(
      snpsift_filter_indels.out.filt_vcfs,
      tx_list_to_indel_norm_and_tumor_cds_fastas.out.somatic_vcfs_w_var_tx_seqs,
      lenstools_filter_expressed_variants.out.somatic_transcripts,
      netmhcpan_alleles,
      pyclonevi_write_results_file.out.pcvi_results,
      procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
      samtools_index.out.bams_and_bais.filter{ it[1] =~ 'ar-' },
      params.manifest_to_lens$filtered_indels_to_class1_neos$gtf,
      params.manifest_to_lens$filtered_indels_to_class1_neos$netmhcpan_parameters,
      params.manifest_to_lens$filtered_indels_to_class1_neos$lenstools_add_indel_metadata_parameters,
      params.manifest_to_lens$filtered_indels_to_class1_neos$netmhcstabpan_parameters,
      params.manifest_to_lens$filtered_indels_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_indels_to_class1_neos$species,
      params.manifest_to_lens$filtered_indels_to_class1_neos$protein_ref)

    // ERV-based antigen prediction
    samtools_coverage_erv(
      samtools_sort_alt.out.bams,
      '')
    if( params.manifest_to_lens$lenstools_filter_expressed_ervs$normal_control_quant ) {
      lenstools_filter_expressed_ervs(
        procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
        params.manifest_to_lens$lenstools_filter_expressed_ervs$normal_control_quant,
        params.manifest_to_lens$lenstools_filter_expressed_ervs_parameters)
      lenstools_filter_ervs_by_rna_coverage(
        lenstools_filter_expressed_ervs.out.expressed_ervs.join(samtools_coverage_erv.out.bam_coverage, by: [0, 1, 2]),
        params.manifest_to_lens$lenstools_filter_ervs_by_rna_coverage_parameters)
    } else if(procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'nr-' }) {
      lenstools_filter_expressed_ervs(
        procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
        procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'nr-' }.map{ it[3] }.first(),
        params.manifest_to_lens$lenstools_filter_expressed_ervs_parameters)
      lenstools_filter_ervs_by_rna_coverage(
        lenstools_filter_expressed_ervs.out.expressed_ervs.join(samtools_coverage_erv.out.bam_coverage, by: [0, 1, 2]),
        params.manifest_to_lens$lenstools_filter_ervs_by_rna_coverage_parameters)
    } else {
      lenstools_filter_expressed_ervs_without_control(
        procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
        params.manifest_to_lens$lenstools_filter_expressed_ervs_without_control$tpm_threshold,
        params.manifest_to_lens$lenstools_filter_expressed_ervs_parameters)
      lenstools_filter_ervs_by_rna_coverage(
        lenstools_filter_expressed_ervs_without_control.out.expressed_ervs.join(samtools_coverage_erv.out.bam_coverage, by: [0, 1, 2]),
        params.manifest_to_lens$lenstools_filter_ervs_by_rna_coverage_parameters)
    } // Keep RAFT happy
    filtered_ervs_to_class1_neos(
      params.dna_ref,
      lenstools_filter_ervs_by_rna_coverage.out.covered_ervs,
      netmhcpan_alleles,
      samtools_index.out.bams_and_bais,
      samtools_sort_alt.out.bams, // This is used by self too, remove _erv suffix.
      procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
      params.manifest_to_lens$filtered_ervs_to_class1_neos$gevequant_general_ref,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$samtools_index_parameters,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$lenstools_get_expressed_ervs_bed_parameters,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$lenstools_make_erv_peptides_parameters,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$netmhcpan_parameters,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$netmhcstabpan_parameters,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$lenstools_add_erv_metadata_parameters,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_ervs_to_class1_neos$species)


    // Viral-based antigen prediction
    lenstools_filter_expressed_viruses(
      unaligned_fqs_to_virdetect_cds_counts.out.viral_cds_counts,
      params.manifest_to_lens$viral_cds_ref,
      params.manifest_to_lens$lenstools_filter_virdetect_by_counts_parameters)
    samtools_sort_viral(
      unaligned_fqs_to_virdetect_cds_counts.out.viral_cds_alns,
      '')
    samtools_coverage_viral(
      samtools_sort_viral.out.bams,
      '')
    lenstools_filter_viruses_by_rna_coverage(
      lenstools_filter_expressed_viruses.out.expressed_viruses.join(samtools_coverage_viral.out.bam_coverage, by: [0, 1, 2]),
      params.manifest_to_lens$lenstools_filter_viruses_by_rna_coverage_parameters)
    samtools_index_viral(
      unaligned_fqs_to_virdetect_cds_counts.out.viral_cds_alns,
      '')
    filtered_viruses_to_class1_neos(
      lenstools_filter_viruses_by_rna_coverage.out.covered_viruses,
      samtools_index_viral.out.bams_and_bais,
      netmhcpan_alleles,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$viral_cds_gff,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$viral_cds_ref,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$viral_pep_ref,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$lenstools_get_expressed_viral_bed_parameters,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$lenstools_make_viral_peptides_parameters,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$netmhcpan_parameters,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$netmhcstabpan_parameters,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$lenstools_add_viral_metadata,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_viruses_to_class1_neos$species)

    // Fusion-based antigen prediction
    filtered_fusions_to_class1_neos(
     manifest_to_rna_procd_fqs.out.procd_fqs.filter{ it[1] =~ 'ar-' },
      procd_fqs_to_starfusion_fusions.out.coding_effect_fusions,
      procd_fqs_to_starfusion_fusions.out.full_fusions,
      netmhcpan_alleles,
      create_phased_germline_tx_vcfs.out.phased_vcfs,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$gtf,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$dna_ref,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$seqtk_subseq_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$netmhcpan_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$netmhcstabpan_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$lenstools_get_fusion_read_count_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$lenstools_add_fusion_metadata_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$bcftools_index_phased_germline_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$lenstools_get_fusion_transcripts_bed_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$samtools_faidx_fetch_parameters,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_fusions_to_class1_neos$species)

    //CTA/Self-antigen prediction
    lenstools_filter_expressed_self_genes(
      procd_fqs_to_star_alns_salmon_counts.out.quants.filter{ it[1] =~ 'ar-' },
      params.manifest_to_lens$gtf,
      params.manifest_to_lens$cta_self_gene_list,
      params.manifest_to_lens$lenstools_filter_expressed_self_parameters)
    filtered_selfs_to_class1_neos(
      samtools_sort_alt.out.bams,
      lenstools_filter_expressed_self_genes.out.expressed_selfs,
      netmhcpan_alleles,
      procd_fqs_to_star_alns_salmon_counts.out.quants,
      snpeff_ann_germline.out.annot_vcfs,
      snpeff_ann.out.annot_vcfs,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$gtf,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$dna_ref,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$cta_self_gene_list,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$samtools_index_parameters,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$lenstools_get_expressed_selfs_bed_parameters,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$samtools_faidx_fetch_parameters,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$bcftools_index_parameters,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$netmhcpan_parameters,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$netmhcstabpan_parameters,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_selfs_to_class1_neos$species)

    //Splice-based antigen prediction
    neosplice_filter_and_summarize_variants(
      bams_to_neosplice.out.neoantigen_results,
      params.manifest_to_lens$peptidome_ref)
    neosplice_filter_and_summarize_variants.out.splice_summaries
      .map{ [it[0], it[2], it[3], it[1], it[4]] }
      .join(samtools_index.out.bams_and_bais, by: [0, 1, 2])
      .map{ [it[0], it[3], it[1], it[2], it[4], it[5], it[6]] }
      .set{ neosplice_summaries_and_bams_and_bais }
    lenstools_get_splice_peptide_read_count(
      neosplice_summaries_and_bams_and_bais,
      params.samps_out_dir)
    neosplice_filter_and_summarize_variants.out.splice_pep_fastas.map{ [it[0], it[2], it[3], it[4]]  }
      .join(netmhcpan_alleles, by: [0, 2])
      .map{ [it[0], it[2], it[1], it[3], it[5]] }
      .set{ splice_peptides_and_alleles }
    netmhcstabpan_rna(
      splice_peptides_and_alleles,
      '')
    antigen_garnish_foreignness_rna(
      lenstools_get_splice_peptide_read_count.out.neosplice_peptide_counts.map{ [it[0], it[2], it[3], it[4]] },
      params.manifest_to_lens$filtered_splice_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_splice_to_class1_neos$species)
    antigen_garnish_dissimilarity_rna(
      lenstools_get_splice_peptide_read_count.out.neosplice_peptide_counts.map{ [it[0], it[2], it[3], it[4]] },
      params.manifest_to_lens$filtered_splice_to_class1_neos$antigen_garnish_data_dir,
      params.manifest_to_lens$filtered_splice_to_class1_neos$species)
    lenstools_get_splice_peptide_read_count.out.neosplice_peptide_counts.map{ [it[0], it[2], it[3], it[4]] }
      .join(antigen_garnish_foreignness_rna.out.foreignness_files, by: [0, 1, 2])
      .join(antigen_garnish_dissimilarity_rna.out.dissimilarity_files, by: [0, 1, 2])
      .join(netmhcstabpan_rna.out.c1_stabs, by: [0, 1, 2])
      .set{ all_splice_metadata_files }

    lenstools_add_splice_metadata(
      all_splice_metadata_files,
      '',
      params.lens_out_dir)

    // Adding TCGA data to snvs and indels
    lenstools_add_tcga_data_snvs(
      filtered_snvs_to_class1_neos.out.metadatas,
      params.manifest_to_lens$tumor_type,
      params.manifest_to_lens$tcga_tx_summary,
      'snv',
      '',
      params.lens_out_dir)

    lenstools_add_tcga_data_indels(
      filtered_indels_to_class1_neos.out.metadatas,
      params.manifest_to_lens$tumor_type,
      params.manifest_to_lens$tcga_tx_summary,
      'indel',
      '',
      params.lens_out_dir)

    filtered_snvs_to_class1_neos.out.metadatas.map{ [it[0], it[4], it[5]] }
      .join(filtered_indels_to_class1_neos.out.metadatas.map{ [it[0], it[4], it[5]] }, by:[0, 1], remainder: true)
      .join(filtered_ervs_to_class1_neos.out.metadatas, by:[0, 1], remainder: true)
      .join(filtered_viruses_to_class1_neos.out.metadatas, by:[0, 1], remainder: true)
      .join(filtered_selfs_to_class1_neos.out.metadatas, by:[0, 1], remainder: true)
      .join(filtered_fusions_to_class1_neos.out.metadatas.map{ [it[0], it[3], it[4]] }, by:[0, 1], remainder: true)
      .join(lenstools_add_splice_metadata.out.metadatas, by:[0, 1], remainder: true)
      .map{ [*it.asList().minus(null)] }
      .map{ [it[0], it[1], it[2..-1]] }
      .set{ joint_reports }

    //Final report
    lenstools_make_lens_report(
      joint_reports,
      '',
      params.lens_out_dir)


    lenstools_make_lens_bed(
        lenstools_make_lens_report.out.reports)

    combine_rel_patient_samples.out.combined_rel_set
      .join(samtools_index.out.bams_and_bais
        .filter{ it[1] =~ 'ar-' }
        .map{ [it[0], it[2], it[1], it[3], it[4]]}, by: [0, 1])
      .set{ all_pat_bams }

    all_pat_bams
      .join(lenstools_make_lens_bed.out.lens_bed, by: [0, 1])
      .set{ all_pat_bams_with_beds }

    bam_subsetter(
      all_pat_bams_with_beds)

    igv_snapshot_automator(
      bam_subsetter.out.subsetted_bams_w_bed)

    // Post-LENS MultiQC

    // Somatic VCF QC
    bcftools_stats_somatic(
      all_raw_somatic_vcfs.mix(all_somatic_vcfs),
      params.manifest_to_lens$bcftools_stats_somatic$bcftools_stats_parameters)

    // Germline VCF QC
    bcftools_stats_germline(
      procd_alns_to_deepvariant_germline_vcfs.out.germline_vcfs.concat(bcftools_filter_germline.out.filtd_vcfs),
      params.manifest_to_lens$bcftools_stats_germline$bcftools_stats_parameters)

    // Alignment QC
    samtools_stats(
      procd_fqs_to_star_alns_salmon_counts.out.alns.concat(procd_fqs_to_procd_alns.out.filt_bams),
      params.manifest_to_lens$samtools_stats$samtools_stats_parameters)

//    multiqc(
//      vcfs_to_vcf_metrics.out.detail_metrics.map{ it -> it[4] }
//        .concat(procd_fqs_to_star_alns_salmon_counts.out.star_logs.map{it -> it[3] })
//        .concat(procd_fqs_to_procd_alns.out.marked_dup_metrics.map{ it -> it[3] })
//        .concat(samtools_stats.out.bam_stats.map{ it -> it[3] })
//        .concat(bcftools_stats_somatic.out.vcf_stats.map{ it -> it[4] })
//        .concat(bcftools_stats_germline.out.vcf_stats.map{ it -> it[4] })
//        .concat(manifest_to_dna_procd_fqs.out.fastqc_zips.map{ it -> [it[3], it[4]] })
//        .concat(manifest_to_rna_procd_fqs.out.fastqc_zips.map{ it -> [it[3], it[4]] })
//        .concat(create_phased_tumor_tx_vcfs.out.phased_stats.map{ it -> it[4] })
//        .concat(create_phased_germline_tx_vcfs.out.phased_stats.map{ it -> it[4] })
//        .collect(),
//     '')
}
