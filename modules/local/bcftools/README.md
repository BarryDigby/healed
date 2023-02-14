# bcftools

BCFtools Nextflow Module

Website: http://samtools.github.io/bcftools/bcftools.html

## Defaults

| workflow               | container                          | cpus              | memory                          |
| ---------------------- | ---------------------------------- | ----------------- | ------------------------------- |
| bcftools_consensus     | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_filter        | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_index         | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_index_somatic | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_isec          | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_merge         | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_norm          | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_query         | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_stats         | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |
| bcftools_stats_somatic | docker://lifebitai/bcftools:1.10.2 | 2 \* task.attempt | 6.GB.plus(4.GB \* task.attempt) |

## Workflows

# bcftools_consensus

Create consensus sequence by applying VCF variants to a reference fasta file.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path(VCF) - VCF File
  path fa - Reference FASTA File
  val parstr - Additional parameters

output:
  tuple => emit: consensus_fastas
    val(pat_name) -  Patient Name
    val(norm_prefix) -  Normal Prefix
    val(tumor_prefix) -  Tumor Prefix
    val(dataset) - Dataset
    path("*.consensus.fa") - Output consensus FASTAs

command:

bcftools consensus ${parstr} -f ${fa} ${vcf} > ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.consensus.fa"
```

# bcftools_filter

Apply fixed-threshold filters.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path(VCF) - VCF File
  path fa - Reference FASTA File
  val parstr - Additional parameters

output:
  tuple => emit: filtd_vcfs
    val(pat_name) -  Patient Name
    val(norm_prefix) -  Normal Prefix
    val(tumor_prefix) -  Tumor Prefix
    val(dataset) - Dataset
    path("*.bfilt.vcf") - Filtered VCFs

command:
bcftools filter ${parstr} -o \${ORIG%.vcf*}.bfilt.vcf ${vcf}
```

# bcftools_index

Creates index for bgzip compressed VCF/BCF files for random access.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
  val parstr - Additional parameters

output:
  tuple => emit: vcfs_w_csis
    val(pat_name) -  Patient Name
    val(prefix) -  Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
    path("*.gz.*i") - VCF Index File

command:
bcftools index ${vcf} ${parstr}
```

# bcftools_index_somatic

Creates index for bgzip compressed VCF/BCF files for random access.
Intended for uses with pair tagged channels (e.g. have both norm_prefix and
tumor_prefix elements).

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal prefix
    val(tumor_prefix) - Tumor prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
  val parstr - Additional parameters

output:
  tuple => emit: vcfs_w_csis
    val(pat_name) -  Patient Name
    val(norm_prefix) -  Normal prefix
    val(tumor_prefix) -  Tumor prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
    path("*.gz.*i") - VCF Index File

command:
bcftools index ${vcf} ${parstr}
```

# bcftools_isec

Creates intersections, unions and complements of VCF files.

```
input:
  tuple => emit: vcfs_w_csis
    val(pat_name) -  Patient Name
    val(norm_prefix) -  Normal prefix
    val(tumor_prefix) -  Tumor prefix
    val(dataset) - Dataset
    path(vcfs) - VCF Files
    path(csis) - VCF Index Files
  val parstr - Additional parameters

output:
  tuple => emit: isec_vcfs
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path("*vcf") - Intersected VCF File

command:
VCFS=\$(ls *vcf.gz); bcftools isec \${VCFS} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.isec.vcf.tmp ${parstr}; echo '##fileformat=VCFv4.2' > ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.isec.vcf; cut -f 1-4 ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.isec.vcf.tmp | sed 's/\\t/\\t\\.\\t/2' >> ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.isec.vcf
```

# bcftools_merge

Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
    path(csi) - VCF Index File
  val parstr - Additional parameters
  val suffix - Output File Suffix

output:
  tuple => emit: merged_vcfs
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Datasets
    val("*vcf") - Merged VCF File

command:
bcftools merge \${VCFS} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.merged.${suffix}.vcf ${parstr}
```

# bcftools_norm

Left-align and normalize indels, check if REF alleles match the reference,
split multiallelic sites into multiple rows; recover multiallelics from
multiple rows.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
  tuple
    path(fa) - Reference FASTA File
    path(faidx) - Reference FASTA Index Files
  val parstr - Additional Parameters

output:
  tuple => emit: norm_vcfs
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path("*norm.vcf") - Normalized VCFa

command:
bcftools norm ${parstr} -f ${fa} ${vcf} > \${BNAME%.vcf*}.norm.vcf
```

# bcftools_query

Extracts fields from VCF or BCF files and outputs them in user-defined format.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
 val parstr - Additional Parameters

output:
  tuple => emit: quots
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
   path("*.qout") - Query Output

command:
bcftools query ${parstr} ${vcf} > ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.qout
```

# bcftools_stats_somatic

Parses VCF or BCF and produces text file stats which is suitable for machine processing.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
 val parstr - Additional Parameters

output:
  tuple => emit: vcf_stats
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    vcf("*vcf_stats") - VCF Statistics File

command:
bcftools stats ${parstr} ${vcf} > \${BNAME%.vcf*}.vcf_stats
```

# bcftools_stats

Parses VCF or BCF and produces text file stats which is suitable for machine processing.

```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(vcf) - VCF File
 val parstr - Additional Parameters

output:
  tuple => emit: vcf_stats
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    vcf("*vcf_stats") - VCF Statistics File

command:
bcftools stats ${parstr} ${vcf} > \${BNAME%.vcf*}.vcf_stats
```
