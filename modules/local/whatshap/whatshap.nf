process whatshap {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label "whatshap_container"
  label "whatshap"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(csi), path(bam), path(bai)
  tuple path(ref), path(idx)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*phased.vcf"), emit: phased_vcfs

  script:
  """
  AVCF=`echo ${vcf}`
  /opt/conda/envs/whp/bin/whatshap phase \
  -o \${AVCF%.vcf*}.phased.vcf  \
  --reference=${ref} \
  ${vcf} \
  ${bam} \ 
  ${parstr}
  """
}


process whatshap_tumor {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label "whatshap_container"
  label "whatshap"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(vcf), path(csi), path(dna_bam), path(dna_bai), path(rna_bam), path(rna_bai)
  tuple path(ref), path(idx)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*phased.vcf"), emit: phased_vcfs

  script:
  """
  AVCF=`echo ${vcf}`
  /opt/conda/envs/whp/bin/whatshap phase \
  --output \${AVCF%.vcf*}.phased.vcf \
  --reference=${ref} \
  ${vcf} \
  ${dna_bam} \
  ${rna_bam} \
  ${parstr} 
  """
}

process whatshap_stats {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label "whatshap_container"
  label "whatshap_stats"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.phased_stats"), emit: phased_stats

  script:
  """
  AVCF=`echo ${vcf}`
  /opt/conda/envs/whp/bin/whatshap stats ${vcf} --tsv=\${AVCF%.phased.vcf}.phased_stats ${parstr}
  """
}
