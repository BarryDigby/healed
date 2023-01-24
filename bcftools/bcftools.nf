#!usr/bin/env nextflow

include { htslib_bgzip } from '../htslib/htslib.nf'
include { samtools_faidx_fetch } from '../samtools/samtools.nf'

process bcftools_sort {
// Creates index for bgzip compressed VCF/BCF files for random access.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: vcfs_w_csis
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path("*sorted*") - Output Sorted VCF

// require:
//   VCFS
//   params.bcftools$bcftool_sort_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_sort'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_sort"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*sorted*"), emit: sorted_vcfs

  script:
  """
  export TMPDIR=\${PWD}
  AVCF=\$(echo ${vcf})
  AVCF_PREFIX=\$(echo \${AVCF%.gz} | rev | cut -f 2 -d '.' | rev)
  AVCF_SUFFIX=\$(echo \${AVCF%.gz} | rev | cut -f 1 -d '.' | rev)
  NEW_OUT="\${AVCF_PREFIX}.sorted.\${AVCF_SUFFIX}"
  bcftools sort ${vcf} ${parstr} > \${NEW_OUT}
  """
}


process bcftools_consensus {
// Create consensus sequence by applying VCF variants to a reference fasta file.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//     path(csi) - CSI File
//     path(fa) - Reference FASTA
//   val parstr - Additional parameters
//
// output:
//   tuple => emit: consensus_fastas
//     val(pat_name) -  Patient Name
//     val(norm_run) -  Normal Run Name
//     val(tumor_run) -  Tumor Run Name
//     val(dataset) - Dataset
//     path("*.consensus.fa") - Output consensus FASTAs

// require:
//   VCFS
//   params.bcftools$bcftools_consensus$dna_ref
//   params.bcftools$bcftools_consensus_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_consensus'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_consensus"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(csi), path(fa)
  val suffix
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*consensus.fa"), emit: consensus_fastas

  script:
  """
  bcftools consensus ${parstr} -f ${fa} ${vcf} > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.${suffix}.consensus.fa
  """
}

process bcftools_filter {
// Apply fixed-threshold filters.
// Thresholds should be passed as string to parstr channel.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(VCF) - VCF
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: filtd_vcfs
//     val(pat_name) -  Patient Name
//     val(norm_run) -  Normal Run Name
//     val(tumor_run) -  Tumor Run Name
//     val(dataset) - Dataset
//     path("*.bfilt.vcf") - Output Filtered VCFs

// require:
//   VCFS
//   params.bcftool$bcftool_filter_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_filter'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_filter"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.bfilt.vcf"), emit: filtd_vcfs

  script:
  """
  ORIG="${vcf}"
  bcftools filter ${parstr} -o \${ORIG%.vcf*}.bfilt.vcf ${vcf}
  """
}


process bcftools_index {
// Creates index for bgzip compressed VCF/BCF files for random access.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: vcfs_w_csis
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//     path("*.*i") - VCF Index

// require:
//   VCFS
//   params.bcftool$bcftool_index_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_index'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_index"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path(vcf), path("*.*i"), emit: vcfs_w_csis

  script:
  """
  bcftools index ${vcf} ${parstr}
  """
}


process bcftools_index2 {
// Creates index for bgzip compressed VCF/BCF files for random access.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//   val parstr - Additional parameters
//
// output:
//   tuple => emit: vcfs_w_csis
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//     path("*.*i") - VCF Index

// require:
//   VCFS
//   params.bcftool$bcftool_index_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_index'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_index"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path(vcf), path("*.*i"), emit: vcfs_w_csis

  script:
  """
  bcftools index ${vcf} ${parstr}
  """
}


process bcftools_index_somatic {
// Creates index for bgzip compressed VCF/BCF files for random access.
// Intended for uses with pair tagged channels (e.g. have both norm_run and
// tumor_run elements).
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: vcfs_w_csis
//     val(pat_name) -  Patient Name
//     val(norm_run) -  Norm Run Name
//     val(tumor_run) -  Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//     path("*.gz.*i") - VCF Index

// require:
//   VCFS
//   params.bcftool$bcftool_index_somatic_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_index'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_index"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path("*.gz.*i"), emit: vcfs_w_csis

  script:
  """
  bcftools index ${vcf} ${parstr}
  """
}


process bcftools_isec {
// Creates intersections, unions and complements of VCF files.
//
// input:
//   tuple => emit: vcfs_w_csis
//     val(pat_name) -  Patient Name
//     val(norm_run) -  Normal Run Name
//     val(tumor_run) -  Tumor Run Name
//     val(dataset) - Dataset
//     path(vcfs) - VCFs
//     path(csis) - CSI Files
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: isec_vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("*vcf") - Intersected VCFs

// require:
//   vcfs_w_csis
//   params.bcftools$bcftools_isec_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_isec'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_isec"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcfs), path(csis)
  val(parstr)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*vcf'), emit: isec_vcfs

  script:
  """
  VCFS=\$(ls *vcf.gz)

  bcftools isec \${VCFS} -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.isec.vcf.tmp ${parstr}
  echo '##fileformat=VCFv4.2' > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.isec.vcf
  cut -f 1-4 ${dataset}-${pat_name}-${norm_run}_${tumor_run}.isec.vcf.tmp | sed 's/\\t/\\t\\.\\t/2' >> ${dataset}-${pat_name}-${norm_run}_${tumor_run}.isec.vcf
  sed -i '2s;^;#CHROM\\tPOS\\tID\\tREF\\tALT\\\n;' *vcf
  """
}


process bcftools_merge {
// Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//     path(csi) - CSI File
//   val parstr - Additional Parameters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: merged_vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Datasets
//     val("*vcf") - Merged VCF

// require:
//   vcfs
//   params.bcftools$bcftools_merge_parameters
//   params.bcftools$bcftools_merge_suffix

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_merge'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_merge"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcfs), path(csis)
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*vcf'), emit: merged_vcfs

  script:
  """
  VCFS=\$(ls *vcf.gz)

  bcftools merge \${VCFS} -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.merged.${suffix}.vcf ${parstr}
  """
}


process bcftools_norm {
// Left-align and normalize indels, check if REF alleles match the reference,
// split multiallelic sites into multiple rows; recover multiallelics from
// multiple rows.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//   tuple
//     path(fa) - Reference FASTA
//     path(faidx) - Reference FASTA Index
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: normd_vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("*normd.vcf") - Normalized VCF

// require:
//   VCFS
//   FA_W_FAIDX
//   params.bcftools$bcftools_norm_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_norm'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_norm"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  tuple path(fa), path(faidx)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*normd.vcf"), emit: normd_vcfs

  script:
  """
  BNAME="${vcf}"
  bcftools norm ${parstr} -f ${fa} ${vcf} > \${BNAME%.vcf*}.normd.vcf
  """
}


process bcftools_query {
// Extracts fields from VCF or BCF files and outputs them in user-defined format.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//  val(parstr) - Additionall Parameters
//
// output:
//   tuple => emit: quots
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("*.qout") - Query Output

// require:
//   VCFS
//   params.bcftool$bcftools_query$bcftool_query_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_query'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_query"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.qout"), emit: qouts

  script:
  """
  bcftools query ${parstr} ${vcf} > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.qout
  """
}


process bcftools_stats_somatic {
// Parses VCF or BCF and produces text file stats which is suitable for machine processing.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//  val(parstr) - Additional Parameters
//
// output:
//   tuple => emit: vcf_stats
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     vcf("*vcf_stats") - VCF Statistics

// require:
//   VCFS
//   params.bcftools$bcftools_stats_somatic_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_stats_somatic'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_stats"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*vcf_stats"), emit: vcf_stats

  script:
  """
  BNAME="${vcf}"
  bcftools stats ${parstr} ${vcf} > \${BNAME%.vcf*}.vcf_stats
  """
}


process bcftools_stats {
// Parses VCF or BCF and produces text file stats which is suitable for machine processing.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcf) - VCF
//  val(parstr) - Additional Parameters
//
// output:
//   tuple => emit: vcf_stats
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     vcf("*vcf_stats") - VCF Statistics

// require:
//   VCFS
//   params.bcftools$bcftools_stats_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_stats'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/bcftools_stats"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*vcf_stats"), emit: vcf_stats

  script:
  """
  BNAME="${vcf}"
  bcftools stats ${parstr} ${vcf} > \${BNAME%.vcf*}.vcf_stats
  """
}


process bcftools_filter_germline {
// Apply fixed-threshold filters.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(VCF) - VCF
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: filtd_vcfs
//     val(pat_name) -  Patient Name
//     val(run) -  Normal Run Name
//     val(dataset) - Dataset
//     path("*.bfilt.vcf") - Filtered VCFs

// require:
//   VCFS
//   params.bcftool$bcftool_filter_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_filter'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_filter"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.bfilt.vcf"), emit: filtd_vcfs

  script:
  """
  ORIG="${vcf}"
  bcftools filter ${parstr} -o \${ORIG%.vcf*}.bfilt.vcf ${vcf}
  """
}



process bcftools_consensus_germline {
// Create consensus_germline sequence by applying VCF variants to a reference
// fasta file.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(VCF) - VCF
//     path(fa) - Reference FASTA
//   val parstr - Additional parameters
//
// output:
//   tuple => emit: consensus_germline_fastas
//     val(pat_name) -  Patient Name
//     val(norm_run) -  Normal Run Name
//     val(tumor_run) -  Tumor Run Name
//     val(dataset) - Dataset
//     path("*.consensus_germline.fa") - Output FASTA

// require:
//   VCFS
//   params.bcftools$bcftools_consensus_germline$dna_ref
//   params.bcftools$bcftools_consensus_germline_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_consensus'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_consensus"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf), path(csi), path(fa)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*consensus.fa"), emit: consensus_fastas

  script:
  """
  bcftools consensus ${parstr} -f ${fa} ${vcf} > ${dataset}-${pat_name}-${run}.consensus.fa
  """
}


process bcftools_call {
// Simple variant caller
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(pileup) - Pileup File
//     path(csi) - Pileup Index
//     path(bed) - BED File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: consensus_germline_fastas
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path("*.bcftools.vcf") - Output VCF

// require:
//   VCFS
//   params.bcftools$bcftools_consensus_germline$dna_ref
//   params.bcftools$bcftools_call_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_call'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_call"

  input:
  tuple val(pat_name), val(run), val(dataset), path(pileup), path(csi), path(bed)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*vcf*"), emit: vcfs

  script:
  """
  bcftools call --regions-file ${bed} ${parstr} ${pileup} > ${dataset}-${pat_name}-${run}.bcftools.vcf
  """
}


process bcftools_mpileup {
// Generate VCF or BCF containing genotype likelihoods for one or multiple
// alignment (BAM or CRAM) files. This is based on the original samtools
// mpileup command (with the -v or -g options) producing genotype likelihoods
// in VCF or BCF format, but not the textual pileup output.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - BAM File
//     path(bai) - BAI File
//     path(ned) - BED File
//   path fa - Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: consensus_germline_fastas
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path("*.pileup") - Output Pileups

// require:
//   VCFS
//   params.bcftools$bcftools_consensus_germline$dna_ref
//   params.bcftools$bcftools_mpileup_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bcftools_container'
  label 'bcftools_mpileup'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/bcftools_mpileup"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path(bed)
  path fa
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.pileup"), path(bed), emit: pileups

  script:
  """
  bcftools mpileup ${parstr} -R ${bed} --fasta-ref ${fa} ${bam} > ${dataset}-${pat_name}-${run}.pileup
  """
}


process bcftools_concat {
// Merge multiple VCF/BCF files from non-overlapping sample sets to create one
// multi-sample file.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcfs) - VCFs
//     path(csis) - VCF Indices
//   val parstr - Additional Parameters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: concatd_vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Datasets
//     val("*vcf") - Merged VCF

// require:
//   vcfs
//   params.bcftools$bcftools_concat_parameters
//   params.bcftools$bcftools_concat_suffix

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_concat'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_concat"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcfs), path(csis)
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*vcf'), emit: concatd_vcfs

  script:
  """
  VCFS=\$(ls -v *vcf* | grep -v csi)

  bcftools concat \${VCFS} -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.concatd.${suffix}.vcf ${parstr} --threads ${task.cpus}
  """
}


workflow bcftools_simple_germline_consensus {
// Creates a germline consensus sequence using BAMs, a genomic reference, and a
// BED defining regions of interest.
//
// Intended to be used for cases where variants cannot be detected easily using
// traditional variant callers (e.g. DeepVariant).
//
// Workflow:
// Combine BAMs with sample-specific beds
//  `-> bcftools mpileup
//       `-> bcftools sort (on bcftools mpileup output)
//            `> bcftools index (on bcftools sort output)
//                `> bcftools call (call variants on sorted, indexed mpileup output)
//                    `> samtools faidx (to get genomic sequence of BED regions)
//                       `> bcftools index (on bcftools call VCF)
//                          `> bcftools consensus (on targeted FASTAs with VCF)
// take:
//   bams - BAM files
//   ref - Reference FASTA
//   beds - Regions of interet BED files
//
// emit:
//   consensus_fastas - Consensus FASTAs (sequences defined by bounds in BED,
//                                        variants from BAM included)
//   vcfs - Intermediate Variant Call Files emitted by BCFtools call function.

// require:
//   BAMS
//   params.bcftools$bcftools_simple_germline_consensus$dna_ref
//   BEDS
  take:
    bams
    ref
    beds
  main:
    bams
      .join(beds, by: [0, 1, 2])
      .set{ bams_and_beds }
    bcftools_mpileup(
      bams_and_beds,
      ref,
      '-Oz')
    bcftools_sort(
      bcftools_mpileup.out.pileups.map{ [it[0], it[1], it[2], it[3]] }, // Stripping out the BED file.
      '-Oz')
    bcftools_index(
      bcftools_sort.out.sorted_vcfs,
      '')
    bcftools_index.out.vcfs_w_csis
      .join(beds, by: [0, 1, 2])
      .set{ pileups_and_beds }
    bcftools_call(
      pileups_and_beds,
      '-c --variants-only -Oz')
    samtools_faidx_fetch(
      beds,
      ref,
      'expressed_hervs',
      '')
    bcftools_index2(
      bcftools_call.out.vcfs,
      '')
    bcftools_index2.out.vcfs_w_csis
      .join(samtools_faidx_fetch.out.fetched_fastas, by: [0, 1, 2])
      .set{ vcfs_w_csis_w_refs }
    bcftools_consensus_germline(
      vcfs_w_csis_w_refs,
      '-H R')
  emit:
    consensus_fastas = bcftools_consensus_germline.out.consensus_fastas
    vcfs = bcftools_call.out.vcfs
}


process bcftools_view {
// View, subset and filter VCF or BCF files by position and filtering
// expression.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcfs) - VCFs
//     path(csis) - VCF Indices
//   val parstr - Additional Parameters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: viewed_vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Datasets
//     val("*vcf") - Merged VCF

// require:
//   VCFS
//   params.bcftools$bcftools_view_parameters
//   params.bcftools$bcftools_view_suffix

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_view'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_view"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(csi)
  val suffix
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*vcf'), emit: viewed_vcfs

  script:
  """
  AVCF="${vcf}"
  bcftools view --samples ${dataset}-${pat_name}-${tumor_run} ${vcf} > \${AVCF%.vcf.gz}.${suffix}.vcf
  """
}


process bcftools_reheader {
// Modify header of VCF/BCF files, change sample names.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(vcfs) - VCFs
//   val parstr - Additional Parameters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: viewed_vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Datasets
//     val("*reh.vcf.gz") - Reheadered VCF

// require:
//   VCFS
//   params.bcftools$bcftools_reheader_parameters
//   params.bcftools$bcftools_reheader_suffix

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'bcftools_container'
  label 'bcftools_reheader'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bcftools_header"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*reh.vcf.gz'), emit: reheadered_vcfs

  script:
  """
  SAMPLE=`grep CHROM ${vcf} | rev | cut -f 1 | rev`
  echo \${SAMPLE}
  AVCF="${vcf}"

  bcftools reheader --samples <(echo \${SAMPLE} ${dataset}-${pat_name}-${norm_run}_${tumor_run}) ${vcf} > \${AVCF%.vcf*}.reh.vcf.gz
  """
}
