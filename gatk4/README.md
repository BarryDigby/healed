# gatk4

GATK4 Nextflow Module

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| gatk_apply_bqsr | docker://broadinstitute/gatk:4.1.4.1 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| gatk_base_recalibrator | docker://broadinstitute/gatk:4.1.4.1 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| gatk_filter_mutect_calls | docker://broadinstitute/gatk:4.1.4.1 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| gatk_haplotypecaller | docker://broadinstitute/gatk:4.1.4.1 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| gatk_index_feature_file | docker://broadinstitute/gatk:4.1.4.1 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| gatk_mutect | docker://broadinstitute/gatk:4.1.4.1 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |

## Workflows

`gatk_haplotypecaller`

Runs gatk HaplotypeCaller on a sample
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    path(bam) - Alignment File
  tuple
    path(fa) - Reference FASTA
    path(idx_files) - Index Files
    path(dict_file) - Dict File
  params - Additional Parameters
  suffix - Output File Suffix
  out_dir - Output Directory

output:
  tuple => emit: germline_vcfs
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path("*vcf") - Variant Call File
```

`gatk_mutect2_matched`

Runs gatk Mutect2 on paired samples
```
input:
  tuple
    val(pat_name) - Patient Name
    val(dataset) - Dataset
    val(norm_prefix) - Normal FASTQ Prefix
    path(norm_bam) - Normal BAM
    path(norm_bai) - Normal BAI
    val(tumor_prefix) - Tumor FASTQ Prefix
    path(tumor_bam) - Tumor BAM
    path(tumor_bai) - Tumor BAI
  tuple
    path(fa) - Reference FASTA
    path(faidx_file) - FAIDX File
    path(dict_file) - Dict File
  val parstr - Additional Paramters
  val suffix - Output File Suffix
  tuple
    path(pon_vcf) - Panel of Normals VCF
    path(pon_vcf_tbi) - Panel of Normals VCF Index
  tuple
    path(af_vcf) - Allele Frequency VCF
    path(af_vcf_tbi) - Allele Frequency VCF Index

output:
  tuple => emit: vcfs
    val(pat_name) - Patient Name
    val(dataset) - Dataset
    path("*.vcf") - VCF
  tuple => emit: vcfs_w_stats
    val(pat_name) - Patient Name
    val(dataset) - Dataset
    path("*.vcf") - VCF
    path("*stats*") - VCF Statistics
````

`gatk_base_recalibrator`

Runs gatk BaseRecalibrator
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(aln) - Alignment File
  tuple
    path(fa) - Reference FASTA
    path(idx_files) - Index Files
    path(dict_file) - Dict File
  tuple
    path(known_sites) - Known Sites Reference VCF
    path(known_sites_index) - Known Sites Reference VCF Index
  val params - Additional Parameters
  val out_dir - Output Directory

output:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(aln) - Alignment File (same as input)
    path("*grp") - Recalibration Matrix
```

`gatk_apply_bqsr`

Runs gatk ApplyBQSR
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(aln) - Alignment File
    path(grp) - Recalibration Matrix
  tuple
    path(fa) - Reference FASTA
    path(idx1) - ???
    path(idx2) - ???
  val parstr - Additional Parameters

output:
  tuple => emit: alns
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path("*bam") - Alignment Files
```

`gatk_index_feature_file`

Runs gatk IndexFeatureFile
```
input:
  path feature_file - Feature File
  val parstr - Additional Parameters

output:
  tuple => emit: ff_w_index
    path(feature_file) - Feature File
    path("${feature_file}.tbi"}) - Feature File Index
```

`gatk_filter_mutect_calls`

Runs gatk FilterMutectCalls on VCFs.
```
input:
  tuple
    val(pat_name) - Patient Name
    val(dataset) - Dataset
    val(norm_prefix) - Normal Prefix
    path(norm_bam) - Normal BAM
    path(norm_bai) - Normal BAI
    val(tumor_prefix) - Tumor Prefix
    path(tumor_bam) - Tumor BAM
    path(tumor_bai) - Tumor BAI
  tuple
    path(fa) - Reference FASTA
    path(faidx_file) - FAIDX File
    path(dict_file) - Dict File
  val parstr - Additional Paramters
  val suffix - Output File Suffix

output:
  tuple => emit: vcfs
    val(pat_name) - Patient Name
    val(dataset) - Dataset
    path("*.vcf")
```
