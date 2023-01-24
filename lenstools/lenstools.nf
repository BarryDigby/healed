#!/usr/bin/env nextflow

/* SNV Processes */

process lenstools_make_snv_peptides {
// Note: This creates both mutant and reference sequences for the purposes of
// calculating agretopicity. It currently does _not_ account for any of the
// patient's germline variants (but it should).

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_snv_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  path tx_aa
  path tx_cds

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.mt_aa.fa"), emit: mutant_fastas
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.wt_aa.fa"), emit: wildtype_fastas
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.wt_aa.fa"), path("*.mt_aa.fa"), emit: peptide_fastas

  script:
  """
  python /opt/lenstools/lenstools.py make-snv-peptides \
  -v ${vcf} \
  -t ${tx_aa} \
  -c ${tx_cds} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.snvs.mt_aa.fa \
  -w ${dataset}-${pat_name}-${norm_run}_${tumor_run}.snvs.wt_aa.fa
  """
}


process lenstools_make_snv_peptides_context {
// Note: This creates both mutant and reference sequences for the purposes of
// calculating agretopicity.

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_snv_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path(var_tx_seqs), path(expressed_txs)
  path gtf
  path pep_ref

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.mt_aa.fa"), emit: mutant_fastas
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.wt_aa.fa"), emit: wildtype_fastas
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.wt_aa.fa"), path("*.mt_aa.fa"), emit: peptide_fastas
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.mt_nt.fa"), emit: nt_fastas

  script:
  """
  grep -f ${expressed_txs} ${gtf} | grep "	CDS	" | cat > subbed_gtf

  python /opt/lenstools/lenstools.py make-snv-peptides-context \
  -sv ${somatic_vcf} \
  -g subbed_gtf \
  -st ${expressed_txs}\
  -vts ${var_tx_seqs} \
  --pep-ref ${pep_ref} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mt_aa.fa \
  --nt-output ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mt_nt.fa \
  -w ${dataset}-${pat_name}-${norm_run}_${tumor_run}.wt_aa.fa \
  --debug-output ${dataset}-${pat_name}-${norm_run}_${tumor_run}.debug
  """
}


process lenstools_add_snv_metadata {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}/lenstools_add_snv_metadata"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(netmhcpan_input), path(pcvi), path(mutant_fastas), path(quants), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  path gtf
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path("*metadata.txt"), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-snv-metadata \
  -b ${netmhcpan_input} \
  -m ${mutant_fastas} \
  -q ${quants} \
  -g ${gtf} \
  -c ${pcvi} \
  -s ${netmhcstabpan} \
  -f ${ag_foreign} \
  -d ${ag_dissim} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}_${rna_run}.snv.metadata.txt
  """
}

/* InDel Processes */

process lenstools_make_indel_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_indel_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  path tx_aa
  path tx_cds

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.aa.fa"), emit: peptide_fastas

  script:
  """
  python /opt/lenstools/lenstools.py make-indel-peptides \
  -v ${vcf} \
  -t ${tx_aa} \
  -c ${tx_cds} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.indels.aa.fa
  """
}

process lenstools_make_indel_peptides_context {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_indel_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path(var_tx_seqs), path(expressed_txs)
  path gtf
  path pep_ref

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.aa.fa"), emit: peptide_fastas
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.nt.fa"), emit: nt_fastas

  script:
  """
  grep -f ${expressed_txs} ${gtf} | grep "	CDS	" | cat > subbed_gtf

  python /opt/lenstools/lenstools.py make-indel-peptides-context \
  -sv ${somatic_vcf} \
  -g subbed_gtf \
  -st ${expressed_txs}\
  --pep-ref ${pep_ref} \
  -vts ${var_tx_seqs} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.aa.fa \
  --nt-output ${dataset}-${pat_name}-${norm_run}_${tumor_run}.nt.fa \
  --debug-output ${dataset}-${pat_name}-${norm_run}_${tumor_run}.debug
  """
}


process lenstools_add_indel_metadata {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}/lenstools_add_indel_metadata"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(netmhcpan_input), path(pcvi), path(mutant_fasta), path(quants), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  path gtf
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path("*metadata.txt"), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-indel-metadata \
  -b ${netmhcpan_input} \
  -m ${mutant_fasta} \
  -q ${quants} \
  -g ${gtf} \
  -c ${pcvi} \
  -s ${netmhcstabpan} \
  -f ${ag_foreign} \
  -d ${ag_dissim} \
  -o "${dataset}-${pat_name}-${norm_run}_${tumor_run}_${rna_run}.indel.metadata.txt"
  """
}

/* ERV Processes */

process lenstools_filter_expressed_ervs {

  label "lenstools_r_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_filter_expressed_ervs"

  input:
  tuple val(pat_name), val(run), val(dataset), path(quants)
  path norm_control
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.expressed_ervs.txt'), emit: expressed_ervs

  script:
  """
  Rscript /opt/lenstools/lenstools.R \
  ${parstr} \
  -t ${quants} \
  -n ${norm_control} \
  -o ${dataset}-${pat_name}-${run}.expressed_ervs.txt
  """
}

process lenstools_filter_ervs_by_rna_coverage {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_filter_ervs_by_rna_coverage"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_ervs), path(coverage)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.expressed_ervs.rc.txt'), optional: true, emit: covered_ervs

  script:
  """
  python /opt/lenstools/lenstools.py check-erv-rna-coverage \
  ${parstr} \
  -c ${coverage} \
  -e ${expressed_ervs} \
  -o ${dataset}-${pat_name}-${run}.expressed_ervs.rc.txt
  """
}


process lenstools_filter_viruses_by_rna_coverage {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_filter_viruses_by_rna_coverage"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_viruses), path(coverage)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.expressed_viruses.rc.txt'), optional: true, emit: covered_viruses

  script:
  """
  python /opt/lenstools/lenstools.py check-virus-rna-coverage \
  ${parstr} \
  -c ${coverage} \
  -e ${expressed_viruses} \
  -o ${dataset}-${pat_name}-${run}.expressed_viruses.rc.txt
  """
}


process lenstools_get_expressed_ervs_bed {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_expressed_ervs_bed"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_ervs)
  path geve_general
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.expressed_ervs.bed'), emit: expressed_ervs_beds

  script:
  """
  python /opt/lenstools/lenstools.py get-expressed-ervs-bed \
  ${parstr} \
  -e ${expressed_ervs} \
  -r ${geve_general} \
  -o ${dataset}-${pat_name}-${run}.expressed_ervs.bed
  """
}


process lenstools_get_expressed_selfs_bed {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_expressed_selfs_bed"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_selfs)
  path gtf
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.expressed_selfs.bed'), emit: expressed_selfs_beds

  script:
  """
  python /opt/lenstools/lenstools.py get-expressed-selfs-bed \
  ${parstr} \
  -e ${expressed_selfs} \
  -g ${gtf} \
  -o ${dataset}-${pat_name}-${run}.expressed_selfs.bed
  """
}


process lenstools_get_expressed_viral_bed {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_expressed_viral_bed"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_viruses)
  path viral_cds_ref
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.expressed_viruses.bed'), emit: expressed_viral_beds

  script:
  """
  python /opt/lenstools/lenstools.py get-expressed-viral-bed \
  ${parstr} \
  -e ${expressed_viruses} \
  -r ${viral_cds_ref} \
  -o ${dataset}-${pat_name}-${run}.expressed_viruses.bed
  """
}


process lenstools_make_erv_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_make_erv_peptides"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_ervs), path(patient_fasta)
  path geve_general
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.ervs.peptides.fa'), path('*.ervs.nt.fa'), emit: erv_peptides
  tuple val(pat_name), val(run), val(dataset), path('*.ervs.nt.fa'), emit: erv_nts

  script:
  """
  python /opt/lenstools/lenstools.py make-erv-peptides \
  ${parstr} \
  -e ${expressed_ervs} \
  -r ${patient_fasta} \
  -g ${geve_general} \
  -o ${dataset}-${pat_name}-${run}.ervs.peptides.fa \
  -n ${dataset}-${pat_name}-${run}.ervs.nt.fa
  """
}


process lenstools_add_erv_metadata {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_add_erv_metadata"

  input:
  tuple val(pat_name), val(run), val(dataset), path(netmhcpan), path(peptides), path(nts), path(quants), path(vcf), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  path geve_data
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(dataset), path('*.ervs.metadata.txt'), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-erv-metadata ${parstr} \
  -p ${peptides} \
  -b ${netmhcpan} \
  -q ${quants} \
  -n ${nts} \
  -v ${vcf} \
  -d ${geve_data} \
  -s ${netmhcstabpan} \
  -f ${ag_foreign} \
  -i ${ag_dissim} \
  -o ${dataset}-${pat_name}-${run}.ervs.metadata.txt
  """
}

/* CTA/Self-antigen Processes */

process lenstools_filter_expressed_self_genes {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_filter_expressed_self_genes"

  input:
  tuple val(pat_name), val(run), val(dataset), path(quants)
  path gtf
  path gene_list
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.expressed_selfs.txt"), emit: expressed_selfs

  script:
  """
  python /opt/lenstools/lenstools.py filter-expressed-self-genes ${parstr} \
  -g ${gene_list} \
  --gtf ${gtf} \
  -q ${quants} \
  -o ${dataset}-${pat_name}-${run}.expressed_selfs.txt
  """
}


process lenstools_make_self_antigen_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_make_self_antigen_peptides"

  input:
  //Currently running without run...
  tuple val(pat_name), val(run), val(dataset), path(consensus_fastas), path(expressed_selfs)
  path gff

  output:
  tuple val(pat_name), val(run), val(dataset), path("*self.pep.fa"), emit: self_antigen_peptides
  tuple val(pat_name), val(run), val(dataset), path("*self.nt.fa"), emit: self_antigen_nts

  script:
  """
  python /opt/lenstools/lenstools.py make-self-antigen-peptides -s ${consensus_fastas} \
  -e ${expressed_selfs} \
  -g ${gff} \
  -o ${dataset}-${pat_name}.self.pep.fa \
  -n ${dataset}-${pat_name}.self.nt.fa
  """
}


process lenstools_add_self_antigen_metadata {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_add_self_antigen_metadata"

  input:
  tuple val(pat_name), val(run), val(dataset), path(binding_affinities), path(quants), path(fasta), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  path gtf
  path gene_list
  val lens_out_dir

  output:
  tuple val(pat_name), val(dataset), path("*self_antigen.metadata.txt"), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-self-antigen-metadata \
  -q ${quants} \
  -b ${binding_affinities} \
  -g ${gtf} \
  -l ${gene_list} \
  -f ${fasta} \
  -s ${netmhcstabpan} \
  -r ${ag_foreign} \
  -d ${ag_dissim} \
  -o ${dataset}-${pat_name}_${run}.self_antigen.metadata.txt
  """
}

/* Viral Processes */

process lenstools_filter_expressed_viruses {

  label 'lenstools_container'
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_filter_expressed_viruses"

  input:
  tuple val(pat_name), val(run), val(dataset), path(viral_quants)
  path viral_ref
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*expressed_viruses.txt'), optional: true, emit: expressed_viruses

  script:
  """
  python /opt/lenstools/lenstools.py filter-expressed-viruses \
  ${parstr} \
  --viral-quants ${viral_quants} \
  --viral-ref ${viral_ref} \
  -o ${dataset}-${pat_name}-${run}.expressed_viruses.txt
  """
}


process lenstools_get_viral_cds_expression {

  label 'lenstools_container'
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_viral_cds_expression"

  input:
  tuple val(pat_name), val(run), val(dataset), path(expressed_viruses), path(viral_cds_counts)
  path viral_cds_ref
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*expressed_viral_cds.txt'), emit: expressed_viruses

  script:
  """
  python /opt/lenstools/lenstools.py get-viral-cds-expression \
  ${parstr} \
  --viral-cds-counts ${viral_cds_counts} \
  --expressed-viruses ${expressed_viruses} \
  -o ${dataset}-${pat_name}-${run}.expressed_viral_cds.txt
  """
}


process lenstools_make_viral_peptides {

  label 'lenstools_container'
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_make_viral_peptides"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fasta)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*viral.pep.fa'), emit: viral_peptides

  script:
  """
  python /opt/lenstools/lenstools.py make-viral-peptides \
  -f ${fasta} \
  -o ${dataset}-${pat_name}-${run}.viral.pep.fa
  """
}


process lenstools_add_viral_metadata {

  label 'lenstools_container'
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_add_viral_metadata"

  input:
  tuple val(pat_name), val(run), val(dataset), path(binding_affinities), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  path virus_cds_ref
  path virus_pep_ref
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(dataset), path('*viral.metadata.txt'), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-viral-metadata \
  ${parstr} \
  --binding-affinities ${binding_affinities} \
  --viral-cds-ref ${virus_cds_ref} \
  --viral-pep-ref ${virus_pep_ref} \
  -s ${netmhcstabpan} \
  -f ${ag_foreign} \
  -d ${ag_dissim} \
  -o ${dataset}-${pat_name}-${run}.viral.metadata.txt
  """
}

/* Fusion Processes */

process lenstools_make_fusion_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_make_fusion_peptides"

  input:
  tuple val(pat_name), val(run), val(run2), val(dataset), path(fusions), path(exon_fas), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  path gtf

  output:
  tuple val(pat_name), val(run), val(run2), val(dataset), path("*.pep.fa"), emit: fusion_peptides
  tuple val(pat_name), val(run), val(run2), val(dataset), path("*.nt.fa"), emit: fusion_nts

  script:
  """
  python /opt/lenstools/lenstools.py make-fusion-peptides \
  -f ${fusions} \
  -o ${dataset}-${pat_name}-${run}.fusion.pep.fa \
  -n ${dataset}-${pat_name}-${run}.fusion.nt.fa
  """
}


process lenstools_make_fusion_nucs {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_make_fusion_nucs"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fusions)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.nuc.fa"), emit: fusion_nucs

  script:
  """
  python /opt/lenstools/lenstools.py make-fusion-nucs \
  -f ${fusions} \
  -o ${dataset}-${pat_name}-${run}.fusion.nuc.fa
  """
}


process lenstools_add_fusion_metadata {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${rna_run}_${norm_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${rna_run}_${norm_run}/lenstools_add_fusion_metadata"

  input:
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path(binding_affinities), path(fasta), path(fusions), path(netmhcstabpan), path(ag_foreign), path(ag_dissim)
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path("*fusion.metadata.txt"), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-fusion-metadata \
  -b ${binding_affinities} \
  -s ${netmhcstabpan} \
  -f ${fusions} \
  -a ${fasta} \
  -r ${ag_foreign} \
  -d ${ag_dissim} \
  -o "${dataset}-${pat_name}-${rna_run}_${norm_run}.fusion.metadata.txt"
  """
}

/* Misc./Shared Processes */

process lenstools_filter_expressed_variants {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_filter_expressed_variants"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(quants)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.expfilt.vcf"), emit: expressed_vcfs
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.somatic_txs"), emit: somatic_transcripts

  script:
  """
  AVCF=`echo ${vcf}`
  gunzip -c ${vcf} > \${AVCF%.gz}

  python /opt/lenstools/lenstools.py filter-expressed-variants \
  ${parstr} \
  -v \${AVCF%.gz} \
  -q ${quants} \
  -o \${AVCF%.vcf*}.expfilt.vcf \
  -s ${dataset}-${pat_name}-${norm_run}_${tumor_run}.somatic_txs

  sed -i 's/^#\\t/1\\t/g' \${AVCF%.vcf*}.expfilt.vcf
  """
}


process lenstools_filter_isolated_variants {
  //Note: This currently does not account for proximal phased heterozygous or
  //hom/het synonymous germline variants. I suspect accounting for proximal
  //somatic variants may be risky since there's no guanrantee they are
  //contained within the same tumor cells.

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_filter_isolated_variants"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path(germline_vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), emit: isolated_vcfs

  script:
  """
  echo "Pass"
  """
}


process lenstools_calculate_agretopicity {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_calculate_agretopicity"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(wt_nmp), path(mt_nmp)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*agreto.netmhcpan.txt"), emit: agreto_netmhcpans

  script:
  """
  MTNMP=`echo ${mt_nmp}`
  python /opt/lenstools/lenstools.py calculate-agretopicity \
  -w ${wt_nmp} \
  -m ${mt_nmp} \
  -o  \${MTNMP%.netmhcpan.txt}.agreto.netmhcpan.txt
  """
}


process lenstools_filter_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_filter_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(binding_affinities)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*filt_peps*"), emit: filtered_peptides

  script:
  """
  BA_TMP=`echo ${binding_affinities}`

  python /opt/lenstools/lenstools.py filter-peptides \
  -i ${binding_affinities} \
  -o  \${BA_TMP%.netmhcpan.txt}.agreto.netmhcpan.txt
  """
}


process lenstools_make_pyclonevi_inputs {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_pyclonevi_inputs"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(candidate_vcf), path(mutect_vcf), path(sequenza_segments), path(sequenza_solutions)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_input"), optional: true, emit: pcvi_inputs

  script:
  """
  echo "mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\terror_rate\ttumour_content" > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_input

  python /opt/lenstools/lenstools.py make-pyclone-vi-inputs \
  -c ${candidate_vcf} \
  -m ${mutect_vcf} \
  -s ${sequenza_segments} \
  --samp-id ${dataset}-${pat_name}-${norm_run}_${tumor_run} \
  --sequenza-solutions ${sequenza_solutions} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_input

  if [ `wc -l ${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_input | cut -f 1 -d ' '` == 1 ]; then 
    mv ${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_input ${dataset}-${pat_name}-${norm_run}_${tumor_run}.empty.pcvi_input
  fi
  """
}

process lenstools_split_pyclonevi_outputs {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_split_pyclonevi_outputs"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(placeholder1), path(placeholder2)
  path pcvi_results
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_results"), emit: pcvi_results

  script:
  """
  head -1 ${pcvi_results} > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_results

  grep "${dataset}-${pat_name}-${norm_run}_${tumor_run}" ${pcvi_results} >> ${dataset}-${pat_name}-${norm_run}_${tumor_run}.pcvi_results
  """
}

process lenstools_get_snv_genomic_context {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_get_snv_genomic_context"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  path tx_cds

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.nuc.fa"), emit: mutant_nucs

  script:
  """
  python /opt/lenstools/lenstools.py get-snv-genomic-context \
  -v ${vcf} \
  -c ${tx_cds} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.snvs.nuc.fa
  """
}


process lenstools_consolidate_multiqc_stats {

  label "lenstools_container"

  input:
  path multiqc_data

  output:
  path "stats.tsv", emit: consolidated_stats

  script:
  """
  python /opt/lenstools/lenstools.py consolidate-multiqc-stats \
  -d ${multiqc_data} \
  -o stats.tsv
  """
}

process add_rna_normals {

  label 'lenstools_container'

  input:
  path manifest
  val run
  val dataset
  val pat_name

  output:
  path "run_data.tsv", emit: run

  script:
  """
  python /opt/lenstools/lenstools.py add-rna-normals \
  -m ${manifest} \
  -p ${run} \
  -d ${dataset} \
  -n ${pat_name} \
  -o run_data.tsv
  """
}

process lenstools_add_splice_metadata {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_add_splice_metadata"

  input:
  tuple val(pat_name), val(tumor_run), val(dataset), path(splice_summary), path(ag_foreign), path(ag_dissim), path(netmhcstabpan)
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(dataset), path("*metadata.txt"), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-splice-metadata \
  -s ${splice_summary} \
  -f ${ag_foreign} \
  -d ${ag_dissim} \
  -b ${netmhcstabpan} \
  -o ${dataset}-${pat_name}-${tumor_run}.splice.metadata.txt
  """
}

process lenstools_make_lens_report {

  label "lenstools_container"
  tag "${dataset}/${pat_name}"
  publishDir "${params.lens_out_dir}/${dataset}/${pat_name}/lenstools_make_lens_reports"

  input:
  tuple val(pat_name), val(dataset), path(metadata_reports)
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(dataset), path("*lens_report.txt"), emit: reports

  script:
  """
  python /opt/lenstools/lenstools.py make-lens-report \
  -d \${PWD} \
  -o ${dataset}-${pat_name}.lens_report.txt
  """
}


process lenstools_add_tcga_data {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}/lenstools_add_tcga_data"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(metadata_report)
  val tumor_type
  path tcga_tx_summ
  val var_type
  val parstr
  val lens_out_dir

  output:
  tuple val(pat_name), val(dataset), val(norm_run), val(tumor_run), val(rna_run), path("*metadata.tcga.txt"), emit: metadatas

  script:
  """
  python /opt/lenstools/lenstools.py add-tcga-data \
  -r ${metadata_report} \
  -t ${tumor_type} \
  -s ${tcga_tx_summ} \
  -o ${dataset}-${pat_name}.${var_type}.metadata.tcga.txt
  """
}

process lenstools_get_viral_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_viral_peptide_read_count"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path(consensus_fa), path(netmhcpan)
  val lens_out_dir

  output:
  tuple val(pat_name), val(run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), optional: true, emit: netmhcpan_peptide_counts

  script:
  """
  python /opt/lenstools/lenstools.py get-peptide-read-count \
  -n ${netmhcpan} \
  -b ${bam} \
  -c ${consensus_fa} \
  -o ${dataset}-${pat_name}-${run}.netmhcpan_peptide_counts.tsv
  """
}

process lenstools_get_erv_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_erv_peptide_read_count"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path(consensus_fa), path(netmhcpan)
  val lens_out_dir

  output:
  tuple val(pat_name), val(run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), emit: netmhcpan_peptide_counts

  script:
  """
  python /opt/lenstools/lenstools.py get-peptide-read-count \
  -n ${netmhcpan} \
  -b ${bam} \
  -c ${consensus_fa} \
  -o ${dataset}-${pat_name}-${run}.netmhcpan_peptide_counts.tsv
  """
}


process lenstools_get_self_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_self_peptide_read_count"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path(consensus_fa), path(netmhcpan)
  val lens_out_dir

  output:
  tuple val(pat_name), val(run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), optional: true, emit: netmhcpan_peptide_counts

  script:
  """
  python /opt/lenstools/lenstools.py get-peptide-read-count \
  -n ${netmhcpan} \
  -b ${bam} \
  -c ${consensus_fa} \
  -o ${dataset}-${pat_name}-${run}.netmhcpan_peptide_counts.tsv
  """
}

process lenstools_get_fusion_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${rna_run}_${norm_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${rna_run}_${norm_run}/lenstools_get_fusion_peptide_read_count"

  input:
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path(starfusion_fusions), path(fusion_reads), path(nt_fa), path(netmhcpan)
  val lens_out_dir

  output:
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), optional: true, emit: netmhcpan_peptide_counts

  script:
  """
  python /opt/lenstools/lenstools.py get-fusion-peptide-read-count \
  -n ${netmhcpan} \
  -f ${starfusion_fusions} \
  -r ${fusion_reads} -t ${nt_fa} \
  -o ${dataset}-${pat_name}-${rna_run}_${norm_run}.netmhcpan_peptide_counts.tsv
  """
}

process lenstools_get_splice_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_get_splice_peptide_read_count"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(neosplice_summaries), path(bam), path(bai)
  val lens_out_dir

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*neosplice_peptide_counts.tsv"), optional: true, emit: neosplice_peptide_counts

  script:
  """
  python /opt/lenstools/lenstools.py get-splice-peptide-read-count \
  -n ${neosplice_summaries} \
  -b ${bam} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.neosplice_peptide_counts.tsv
  """
}


//process lenstools_get_snv_peptide_read_count {
//
//
//  label "lenstools_container"
//  tag "${dataset}/${pat_name}/${run}"
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_snv_peptide_read_count"
//
//  input:
//  tuple val(pat_name), val(run), val(dataset), path(netmhcpan), path(bam), path(bai)
//  val lens_out_dir
//
//  output:
//  tuple val(pat_name), val(run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), emit: snv_peptide_counts
//
//  script:
//  """
//  python /opt/lenstools/lenstools.py get-snv-peptide-read-count -n ${netmhcpan} -b ${bam} -o ${dataset}-${pat_name}-${run}.netmhcpan_peptide_counts.tsv
//  """
//}

process lenstools_get_snv_peptide_count {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}/lenstools_get_snv_peptide_count"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(netmhcpan), path(rna_bam), path(rna_bai)
  path tx_cds

  output:
  // This should emit the rna_run as well. Only emitting DNA normal and DNA
  // tumor for the purposes of joining with other outputs for making metadata
  //  reports.
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*snv.netmhcpan"), optional: true, emit: snv_peptide_counts

  script:
  """
  AFASTA=`echo ${peptide_fasta}`

  python /opt/lenstools/lenstools.py get-snv-peptide-count \
  -f ${netmhcpan} \
  -b ${rna_bam} \
  -t ${tx_cds} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}_${rna_run}.snv_peptide_counts.tsv
  """
}

process lenstools_filter_mutant_snv_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_filter_mutant_snv_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(wt_nmp), path(mt_nmp), path(mt_fa)
  val max_peptide_length

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*${norm_run}*wt.filt.netmhcpan.txt"), path("*${tumor_run}*mt.filt.netmhcpan.txt"), emit: filt_netmhcpans

  script:
  """
  python /opt/lenstools/lenstools.py filter-mutant-peptides \
  -w ${wt_nmp} \
  -m ${mt_nmp} \
  -mf ${mt_fa} \
  -wo ${dataset}-${pat_name}-${norm_run}_${tumor_run}.wt.filt.netmhcpan.txt \
  -mo ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mt.filt.netmhcpan.txt \
  -l ${max_peptide_length}
  """
}

process lenstools_filter_mutant_indel_peptides {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_filter_mutant_snv_peptides"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(mt_nmp), path(mt_fa)
  val max_peptide_length

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*${tumor_run}*mt.filt.netmhcpan.txt"), emit: filt_netmhcpans

  script:
  """
  python /opt/lenstools/lenstools.py filter-mutant-peptides \
  -m ${mt_nmp} \
  -mf ${mt_fa} \
  -mo ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mt.filt.netmhcpan.txt \
  -l ${max_peptide_length}
  """
}

process lenstools_get_expressed_transcripts_bed {
//Intended for making beds for transcripts containing somatic variants

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_get_expressed_transcripts_bed"

  input:
  //Should include norm, tumor, and rna vals
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(expressed_transcripts)
  path gtf
  val parstr

  output:
//  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*.expressed_transcripts.bed'), emit: expressed_transcripts_beds
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*tx_beds'), emit: expressed_transcripts_beds

  script:
  """
  mkdir -p ${dataset}-${pat_name}-${norm_run}_${tumor_run}.tx_beds
  while read line; do 
    echo \${line}; 
    grep \${line} *gtf | grep '	CDS	' | cut -f 1,4,5 >> ${dataset}-${pat_name}-${norm_run}_${tumor_run}.tx_beds/${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${line//'"'/}.bed; 
  done < ${expressed_transcripts}
  """
}

process lenstools_get_fusion_transcripts_bed {
//Intended for making beds for transcripts containing somatic variants

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/lenstools_get_fusion_transcripts_bed"

  input:
  //Should include norm, tumor, and rna vals
  tuple val(pat_name), val(run), val(dataset), path(expressed_transcripts)
  path gtf
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.fusion_transcripts.bed'), emit: fusion_transcripts_beds

  script:
  """
  python /opt/lenstools/lenstools.py get-expressed-transcripts-bed \
  ${parstr} \
  -t ${expressed_transcripts} \
  -g ${gtf} \
  -o ${dataset}-${pat_name}-${run}.fusion_transcripts.bed
  """
}

process lenstools_make_variant_specific_vcfs {

  label 'bcftools_container'
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_variant_specific_vcfs"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path(norm_phased_vcf), path(tumor_phased_vcf)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path(all_vcfs), emit: somatic_vcf_and_vcf_spec_vcfs

  script:
  """
  mkdir -p var_vcfs
  cut -f 1,2 ${somatic_vcf} | grep -v "^#" > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.voi_coords

  while read line; do
    echo \${line};
    CHROM=`echo \${line} | cut -f 1 -d ' '`;
    POS=`echo \${line} | cut -f 2 -d ' '`;
    echo "CHROM: \${CHROM} POS: \${POS}";
    bcftools query -i "CHROM=\\"\${CHROM}\\" && POS=\${POS}" -f "%CHROM\t%POS\t[%PS]\n" *tumor.phased.vcf >> voi_metadata;
  done < ${dataset}-${pat_name}-${norm_run}_${tumor_run}.voi_coords

  while read line; do
    echo \${line};
    CHROM=`echo \${line} | cut -f 1 -d ' '`;
    POS=`echo \${line} | cut -f 2 -d ' '`;
    PS=`echo \${line} | cut -f 3 -d ' '`;
    echo \${CHROM} \${POS} \${PS};
    LB=`expr \${POS} - 50`;
    UB=`expr \${POS} + 50`;
    echo \${LB} \${UB};
    PS_STR='';
    if [[ \${PS} =~ ^[0-9]+\$ ]]; then
#      PS_STR=" || CHROM=\\"\${CHROM}\\" && POS > \${LB} && POS < \${UB} && GT=\\"het\\" && PS=\${PS}";
      PS_STR=" || GT=\\"het\\" && PS=\${PS}";
      echo "\$PS_STR";
    else
      echo "Unusable PS";
    fi;
#    FILTER_STR="CHROM=\\"\${CHROM}\\" && POS=\${POS} || CHROM=\\"\${CHROM}\\" && POS > \${LB} && POS < \${UB} && GT=\\"AA\\"\${PS_STR}";
    FILTER_STR="CHROM=\\"\${CHROM}\\" && POS=\${POS} || GT=\\"AA\\"\${PS_STR}";
    echo \${FILTER_STR};
    bcftools filter -i "\${FILTER_STR}" ${norm_phased_vcf} > \${CHROM}_\${POS}.normal.vcf;
    bcftools filter -i "\${FILTER_STR}" ${tumor_phased_vcf} > \${CHROM}_\${POS}.tumor.vcf;
  done < voi_metadata

  mkdir -p all_vcfs
  mv *normal.vcf all_vcfs
  mv *tumor.vcf all_vcfs
  """
}

process lenstools_make_variant_specific_tx_seqs {

  label 'bcftools_container'
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_variant_specific_vcfs"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path(all_vcfs), path(txs_fas)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(somatic_vcf), path("all_var_tx_seqs"), optional: true, emit: somatic_vcfs_w_var_tx_seqs

  script:
  """
  # Want to copy to a local path. Seems cp -P isn't working. Need to fix.
  mkdir -p all_vcfs.local

  for i in `ls all_vcfs/*vcf`; do
    echo \${i}
    OUT=`echo \${i} | rev | cut -f 1 -d '/' | rev`
    bgzip -c \${i} > all_vcfs.local/\${OUT}.gz
  done
  cd all_vcfs.local
  for j in `ls *gz`; do
    echo \${j}
    bcftools index \${j}
  done
  cd ..

  for i in `ls ${txs_fas}/*`; do
    echo \${i};
    ENST=`echo \${i} | rev | cut -f 3,4 -d '.' | rev`;
    echo \${ENST};
    for target_var in `grep \${ENST} *.vcf | cut -f 1,2 | sed 's/\\t/_/g'`;
      do echo \${target_var};
      for vcf in `ls all_vcfs.local/\${target_var}*gz`; do
        VCF_OUT=`echo \${vcf} | rev | cut -f 1 -d '/' | rev`;
        bcftools consensus -H 2pIu --mark-del X --mark-ins lc -f \${i} \${vcf} -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${ENST}_\${VCF_OUT%.vcf*}.fa;
      done;
      if cmp --silent ${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${ENST}_\${target_var}*tumor*fa ${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${ENST}_\${target_var}*norm*fa; then
        rm ${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${ENST}_\${target_var}*fa 
      fi;
    done;
  done

  if ls *tumor.fa 1> /dev/null 2>&1; then
    mkdir -p all_var_tx_seqs
    mv *normal.fa all_var_tx_seqs
    mv *tumor.fa all_var_tx_seqs
  fi
  """
}

process lenstools_get_snv_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}/lenstools_get_snv_peptide_read_count"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(bam), path(bai), path(nt_fa), path(netmhcpan)
  path gtf
  val lens_out_dir

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), optional: true, emit: netmhcpan_peptide_counts

  script:
  """
  grep \\> ${nt_fa} | cut -f 5 -d ' ' | sed 's/TRANSCRIPT://g' > txs_of_interest

   while read line; do
    echo \${line};
    grep \${line} ${gtf} | cat >> subbed_gtf;
  done < txs_of_interest

  if [ -f subbed_gtf ]; then
    python /opt/lenstools/lenstools.py get-snv-peptide-read-count \
    -n ${netmhcpan} \
    -b ${bam} \
    -c ${nt_fa} \
    -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}_${rna_run}.netmhcpan_peptide_counts.tsv \
    -g subbed_gtf
  fi
  """
}

process lenstools_get_indel_peptide_read_count {


  label "lenstools_container"
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}_${rna_run}/lenstools_get_indel_peptide_read_count"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(bam), path(bai), path(nt_fa), path(netmhcpan)
  path gtf
  val lens_out_dir

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path("*netmhcpan_peptide_counts.tsv"), optional: true, emit: netmhcpan_peptide_counts

  script:
  """
  grep \\> ${nt_fa} | cut -f 4 -d ' ' | sed 's/TRANSCRIPT://g' > txs_of_interest

   while read line; do
    echo "Grepping \${line} entires from GTF.";
    grep \${line} ${gtf} | cat >> subbed_gtf;
  done < txs_of_interest

  if [ -f subbed_gtf ]; then
    python /opt/lenstools/lenstools.py get-indel-peptide-read-count \
    -n ${netmhcpan} \
    -b ${bam} \
    -c ${nt_fa} \
    -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}_${rna_run}.netmhcpan_peptide_counts.tsv \
    -g subbed_gtf
  fi
  """
}

process lenstools_make_fusion_specific_vcfs {

  label 'bcftools_container'
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_variant_specific_vcfs"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(fusions)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(all_vcfs), emit: somatic_vcf_and_vcf_spec_vcfs

  script:
  """
  mkdir -p fusion_vcfs

  """
}

process lenstools_make_fusion_specific_tx_seqs {

  label 'bcftools_container'
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/lenstools_make_variant_specific_vcfs"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(exon_fas),  path(germline_vcf), path(csi)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*iupac_germline_vars.fa"), emit: var_tx_seqs

  script:
  """
  bcftools consensus \
  -I \
  --mark-del X \
  --mark-ins lc \
  -f ${exon_fas} \
  ${germline_vcf} \
  -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.fusion_exons.iupac_germline_vars.fa
  """
}

process lenstools_make_fusion_peptides_context {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${rna_run}_${norm_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${rna_run}_${norm_run}/lenstools_make_fusion_peptides_context"

  input:
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path(fusions), path(exon_fas)
  path gtf

  output:
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path("*.pep.fa"), optional: true, emit: fusion_peptides
  tuple val(pat_name), val(rna_run), val(norm_run), val(dataset), path("*.nt.fa"), optional: true, emit: fusion_nts

  script:
  """
  for i in `cut -f 18,20 ${fusions} | grep ENS | sed 's/\\t/\\n/'`; do
    echo \${i}
    grep \${i} ${gtf} | cat >> subbed_gtf
    echo \${i} >> expressed_fusion_txs
  done
 
  if [ -f expressed_fusion_txs ]; then
    python /opt/lenstools/lenstools.py make-fusion-peptides-context \
    -f ${fusions} \
    -e ${exon_fas} \
    -g subbed_gtf \
    -t expressed_fusion_txs \
    -o ${dataset}-${pat_name}-${rna_run}_${norm_run}.fusion.pep.fa \
    -n ${dataset}-${pat_name}-${rna_run}_${norm_run}.fusion.nt.fa
  fi
  """
}

process lenstools_filter_expressed_ervs_without_control {

  label "lenstools_container"
  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(quants)
  val tpm_threshold
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path('*.expressed_ervs.txt'), emit: expressed_ervs

  script:
  """
  python /opt/lenstools/lenstools.py filter-expressed-ervs ${parstr} \
    -q ${quants} \
    -z ${tpm_threshold} \
    -o ${dataset}-${pat_name}-${prefix}.expressed_ervs.txt
  """
}

process lenstools_make_lens_bed {

  label "lenstools_container"
  tag "${dataset}/${pat_name}"

  input:
  tuple val(pat_name), val(dataset), path(lens_report)


  output:
  tuple val(pat_name), val(dataset), path("*bed"), emit: lens_bed

  script:
  """
  python /opt/lenstools/lenstools.py make-lens-bed \
    -d ${lens_report} \
    -s ${dataset}-${pat_name} \
#    -o ${dataset}-${pat_name}.tumor_antigen.bed \
  """
}
