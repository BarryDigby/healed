#!/usr/bin/env nextflow

include { extract_chroms_from_gff } from '../ref_utils/ref_utils.nf'

include { samtools_sort } from '../samtools/samtools.nf'
include { samtools_index } from '../samtools/samtools.nf'

include { combine_patient_samples } from '../utilities/utilities.nf'
include { combine_patient_samples as combine_patient_junctions } from '../utilities/utilities.nf'

process neosplice_augmented_splice_graph_build {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'neosplice_container'
  label 'neosplice_augmented_splice_graph_build'

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai), val(chr)
  path fa
  path gff
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*_graph.json"), emit: graph_jsons

  script:
  """
  CHR=`echo ${chr}`
  BAM_PLACEHOLDER=`echo ${tumor_bam}`
  FA_PLACEHOLDER=`echo ${fa}`
  samtools view -bS ${tumor_bam} \${CHR} > \${BAM_PLACEHOLDER%.bam}.\${CHR}.bam
  samtools index \${BAM_PLACEHOLDER%.bam}.\${CHR}.bam
  samtools faidx ${fa} \${CHR} > \${FA_PLACEHOLDER%.fa}.\${CHR}.fa
  python /NeoSplice/augmented_splice_graph.py build ${parstr} \
      --bam \${BAM_PLACEHOLDER%.bam}.\${CHR}.bam \
      --genome \${FA_PLACEHOLDER%.fa}.\${CHR}.fa \
      --gff  ${gff} \
      --seq \${CHR} \
      --cpus 8 \
      --mem-per-cpu 3
  mv output/*json .
  """
}


process neosplice_get_max_kmer_length {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'neosplice_container'
  label 'neosplice_augmented_splice_graph_build'

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_run), val(tumor_run), path("*max_kmer_len.txt"), emit: max_kmer_len

  script:
  """
  python /NeoSplice/get_max_kmer_length.py \
    --tumor_bam ${tumor_bam} \
    --normal_bam ${norm_bam} \
    > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.max_kmer_len.txt
  """
}


process neosplice_convert_bam_to_fasta {

  tag "${dataset}/${pat_name}/${run}"
  label 'neosplice_container'
  label 'neosplice_convert_bam_to_fasta'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*_R1.fasta"), path("*_R2.fasta"), emit: fastas

  script:
  """
  python /NeoSplice/convert_bam_to_fasta.py \
  --bam_file ${bam} \
  --R1_out ${dataset}-${pat_name}-${run}_R1.fasta \
  --R2_out ${dataset}-${pat_name}-${run}_R2.fasta
  """
}


process neosplice_msbwtis {

  tag "${dataset}/${pat_name}/${run}"
  label 'neosplice_container'
  label 'neosplice_msbwtis'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*tmp"), emit: msbwt_tmps

  script:
  """
  mkdir -p ${dataset}-${pat_name}-${run}.msbwt.tmp
  /msbwt-is/msbwtis ${dataset}-${pat_name}-${run}.msbwt.tmp ${fq1} ${fq2}
  """
}


process neosplice_convert_bwt_format {

  tag "${dataset}/${pat_name}/${run}"
  label 'neosplice_container'
  label 'neosplice_convert_bwt_format'

  input:
  tuple val(pat_name), val(run), val(dataset), path(msbwt_tmp)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.msbwt"), emit: msbwts

  script:
  """
  mkdir -p ${dataset}-${pat_name}-${run}.msbwt
  bash /NeoSplice/convert_BWT_format.bash ${msbwt_tmp} ${dataset}-${pat_name}-${run}.msbwt
  """
}


process neosplice_kmer_search_bwt {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'neosplice_container'
  label 'neosplice_kmer_search_bwt'

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bwt), val(tumor_run), path(tumor_bwt), path(read_length)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_run), val(tumor_run), path("*.tumor_kmers.txt"), emit: tumor_kmers

  script:
  """
  max_kmer_read_len=`cat ${read_length}`
  python /NeoSplice/Kmer_search_bwt.py \
    --tumor_bwt ${tumor_bwt} \
    --normal_bwt ${norm_bwt} \
    --processors 1 --max_length \${max_kmer_read_len} \
    --tumor_threshold 20 \
    --normal_threshold 4  \
    --outdir ${dataset}-${pat_name}-${norm_run}_${tumor_run}.kmers/
  cat ${dataset}-${pat_name}-${norm_run}_${tumor_run}.kmers/Tumor_kmers_* >  ${dataset}-${pat_name}-${norm_run}_${tumor_run}.tumor_kmers.txt
  """
}


process neosplice_search_bam {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'neosplice_container'
  label 'neosplice_search_bam'

  input:
  tuple val(pat_name), val(dataset), val(norm_run), val(tumor_run), path(kmers), path(bam), path(bai)
  val parstr

  output:
  tuple val(pat_name), val(tumor_run), val(dataset), path("*tumor_kmer.bam"), emit: kmer_bams

  script:
  """
  python /NeoSplice/search_bam.py \
  --Kmer_file ${kmers} \
  --input_bam_file ${bam} \
  --out_bam_file ${dataset}-${pat_name}-${norm_run}_${tumor_run}.tumor_kmer.bam
  """
}


process neosplice_get_splice_junctions {

  tag "${dataset}/${pat_name}/${run}"
  label 'neosplice_container'
  label 'neosplice_get_splice_junctions'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*junctions.txt"), emit: junctions

  script:
  """
  python /NeoSplice/get_splice_junctions.py --input_bam ${bam} --out_file ${dataset}-${pat_name}-${run}.junctions.txt
  """
}


process neosplice_kmer_graph_inference {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'neosplice_container'
  label 'neosplice_kmer_graph_inference'

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(tumor_bam), path(tumor_bai), path(tumor_kmer_bam), path(tumor_kmer_bai), path(norm_junctions), path(tumor_junctions), path(hla_calls), path(splice_graph)
  path gff
  path reference_fa
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("${dataset}-${pat_name}-*neoantigen_result"), emit: neoantigen_results

  script:
  """
  mkdir -p tmp
  export NETMHCpan=/netMHCpan-4.1-docker
  export TMPDIR=\${PWD}/tmp

  BAM_PLACEHOLDER=`echo ${tumor_bam}`
  KMER_BAM_PLACEHOLDER=`echo ${tumor_kmer_bam}`
  FA_PLACEHOLDER=`echo ${reference_fa}`

  HLA_ALLELES=`cat ${hla_calls}`
  for i in `ls *json`; do
    CHR=`echo \${i} | cut -f 1 -d '_'`;
    samtools view -bS ${tumor_bam} \${CHR} > \${BAM_PLACEHOLDER%.bam}.\${CHR}.bam
    samtools index \${BAM_PLACEHOLDER%.bam}.\${CHR}.bam
    samtools view -bS ${tumor_kmer_bam} \${CHR} > \${KMER_BAM_PLACEHOLDER%.bam}.\${CHR}.bam
    samtools index \${KMER_BAM_PLACEHOLDER%.bam}.\${CHR}.bam
    samtools faidx ${reference_fa} \${CHR} > \${FA_PLACEHOLDER%.fa}.\${CHR}.fa
    python /NeoSplice/kmer_graph_inference.py --sample ${dataset}-${pat_name}-${norm_run}_${tumor_run} \
      --chromosome \${CHR} \
      --bam_file \${BAM_PLACEHOLDER%.bam}.\${CHR}.bam \
      --gff_file ${gff} \
      --genome_fasta \${FA_PLACEHOLDER%.fa}.\${CHR}.fa \
      --kmer_bam \${KMER_BAM_PLACEHOLDER%.bam}.\${CHR}.bam \
      --splice_graph \${i} \
      --tumor_junction_file ${tumor_junctions} \
      --normal_junction_file ${norm_junctions} \
      --transcript_min_coverage 100 --HLA_I \${HLA_ALLELES} \
      --HLA_II None \
      --netMHCpan_path /netMHCpan-4.1-docker/netMHCpan \
      --netMHCIIpan_path /netMHCIIpan-3.2-docker/netMHCIIpan \
      --outdir ${dataset}-${pat_name}-\${CHR}.neoantigen_result&
  done
  wait

  rm -rf \${PWD}/tmp
  """
}


workflow neosplice {
  take:
    manifest
    bams_bais
    fa
    gff
    hla_calls
    neosplice_augmented_splice_graph_build_parameters
    neosplice_get_max_kmer_length_parameters
    neosplice_convert_bams_to_fasta_parameters
    neosplice_get_splice_junctions_parameters
    neosplice_msbwt_is_parameters
    neosplice_convert_bwt_format_parameters
    neosplice_kmer_search_bwt_parameters
    neosplice_search_bam_parameters
    samtools_sort_parameters
    samtools_index_parameters
    neosplice_kmer_graph_inference_parameters
  main:
    extract_chroms_from_gff(
      gff)
    bams_bais.filter{ it[1] =~ 'ar-' }.set{tumor_bam_bais}
    bams_bais.filter{ it[1] =~ 'nr-' }.set{norm_bam_bais}
    tumor_bam_bais
      .combine(norm_bam_bais)
      .map{ [it[0], it[2], it[6], it[8], it[9], it[1], it[3], it[4]] }
      .set{ rna_bams }
    rna_bams
      .combine(extract_chroms_from_gff.out.chroms_list.splitText())
      .set{ rna_bams_w_chr }
    neosplice_augmented_splice_graph_build(
      rna_bams_w_chr,
      fa,
      gff,
      neosplice_augmented_splice_graph_build_parameters)
    neosplice_get_max_kmer_length(
      rna_bams,
      neosplice_get_max_kmer_length_parameters)
    rna_bams
      .flatMap{ it -> [[it[0], it[2], it[1], it[3], it[4]], [it[0], it[5], it[1], it[6], it[7]]] }
      .set{ sample_level_bams }
    neosplice_convert_bam_to_fasta(
      sample_level_bams,
      neosplice_convert_bams_to_fasta_parameters)
    neosplice_get_splice_junctions(
      sample_level_bams,
      neosplice_get_splice_junctions_parameters)
    neosplice_msbwtis(
      neosplice_convert_bam_to_fasta.out.fastas,
      neosplice_msbwt_is_parameters)
    neosplice_convert_bwt_format(
      neosplice_msbwtis.out.msbwt_tmps,
      neosplice_convert_bwt_format_parameters)
    manifest.filter{ it[5] =~ 'TRUE' }.set{ norm_set }
    manifest.filter{ it[5] =~ 'FALSE' }.set{ tumor_set }
    neosplice_convert_bwt_format.out.msbwts
      .filter{ it[1] =~ 'nr-' }
      .map{ [it[0], it[2], it[1], it[3]] }
      .set{ norm_msbwts }
    neosplice_convert_bwt_format.out.msbwts
      .filter{ it[1] =~ 'ar-' }
      .map{ [it[0], it[2], it[1], it[3]] }
      .set{ tumor_msbwts }
    norm_msbwts
      .join(tumor_msbwts, by: [0, 1])
      .set{ norm_tumor_msbwts }
    norm_tumor_msbwts
      .map{ [it[0], it[1], it[2], it[4], it[3], it[5]] }
      .join(neosplice_get_max_kmer_length.out.max_kmer_len, by: [0, 1, 2, 3])
      .map{ [it[0], it[1], it[2], it[4], it[3], it[5], it[6]] }
      .set{norm_tumor_msbwts_read_len }
    neosplice_kmer_search_bwt(
      norm_tumor_msbwts_read_len,
      neosplice_kmer_search_bwt_parameters)
    rna_bams
      .map{ [it[0], it[1], it[2], it[5], it[6], it[7]] }
      .set{ search_bams_input_component }
    neosplice_kmer_search_bwt.out.tumor_kmers
      .join(search_bams_input_component, by: [0, 1, 2, 3])
      .set{ search_bam_inputs }
    neosplice_search_bam(
      search_bam_inputs,
      neosplice_search_bam_parameters)
    samtools_sort(
      neosplice_search_bam.out.kmer_bams,
      samtools_sort_parameters)
    samtools_index(
      samtools_sort.out.bams,
      samtools_index_parameters)
    bams_bais
      .join(samtools_index.out.bams_and_bais, by: [0, 1, 2])
      .set{ bams_and_kmer_bams }
    neosplice_get_splice_junctions.out.junctions
      .filter{ it[1] =~ 'nr-' }
      .map{ [it[0], it[2], it[1], it[3]] }
      .set{ norm_junctions }
    neosplice_get_splice_junctions.out.junctions
      .filter{ it[1] =~ 'ar-' }
      .map{ [it[0], it[2], it[1], it[3]] }
      .set{ tumor_junctions }
    norm_junctions.join(tumor_junctions, by: [0, 1])
      .set{ norm_tumor_junctions }
    bams_and_kmer_bams
      .join(norm_tumor_junctions, by: 0)
      .map{ [it[0], it[8], it[1], it[2], it[3], it[4], it[5], it[6], it[9], it[11]] }
      .set{ bams_and_kmer_bams_and_junctions }
    bams_and_kmer_bams_and_junctions
      .join(hla_calls, by: 0)
      .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8], it[9], it[12]] }
      .set{ bams_and_kmer_bams_and_junctions_and_hla_alleles }
    bams_and_kmer_bams_and_junctions_and_hla_alleles
      .combine(neosplice_augmented_splice_graph_build.out.graph_jsons, by: [0, 1, 2, 3])
      .set{ bams_and_kmer_bams_and_junctions_and_hla_alleles_and_splice_graphs }
    neosplice_kmer_graph_inference(
      bams_and_kmer_bams_and_junctions_and_hla_alleles_and_splice_graphs,
      gff,
      fa,
      neosplice_kmer_graph_inference_parameters)
    if( params.species =~ /[Hh]uman|hs|HS|[Hh]omo/ ) {
        chrom_count = 25
    } else if( params.species =~ /[Mm]ouse|mm|MM|[Mm]us/ ) {
        chrom_count = 21
    }
    neoantigen_results = neosplice_kmer_graph_inference.out.neoantigen_results
      .groupTuple(by: [0,1,2,3], size:chrom_count)
      .set{ neo_results_by_pat}

  emit:
    neoantigen_results = neo_results_by_pat
}

process neosplice_filter_and_summarize_variants {

  label "neosplice_container"
  label "neosplice_sv_summarization"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/neosplice_filter_and_summarize_variants"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(outcomes_peptides)
  path peptidome_ref

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("SV*txt"), emit: splice_summaries
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*splice.pep.fa"), emit: splice_pep_fastas

  script:
  """
  mkdir -p outcome_peptides
  cp *result/*/*/* outcome_peptides
  python /NeoSplice/splice_variant_neoantigen_summarization_cs.py \
    --ref-dir ${peptidome_ref}/ \
    --outcome-peptides-dir outcome_peptides

  for i in `cut -f 1 SV*txt | tail -n +1`; do
    echo ">\${i}" >> ${dataset}-${pat_name}-${norm_run}_${tumor_run}.splice.pep.fa
    echo "\${i}" >> ${dataset}-${pat_name}-${norm_run}_${tumor_run}.splice.pep.fa
  done
  """
}
