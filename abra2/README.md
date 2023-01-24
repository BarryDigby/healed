# abra2

Abra2 Nextflow Module for use with RAFT

Website: https://github.com/mozack/abra2

Manuscript: https://doi.org/10.1093/bioinformatics/btz033

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| abra2 | docker://aphoid/abra2:2.20 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| abra2_cadabra | docker://aphoid/abra2:2.20 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| abra2_rna | docker://aphoid/abra2:2.20 | 4 * task.attempt |  20.GB.plus(12.GB * task.attempt) |

## Workflows

# abra2

Runs abra2
ABRA is a realigner for next generation sequencing data. It uses localized
assembly and global realignment to align reads more accurately, thus
improving downstream analysis (detection of indels and complex variants in
particular).
```
input:
  tuple
    val(pat_name) - Patient Name
    val(datset) - Dataset
    val(norm_prefix) - Normal Prefix
    val(norm_bam) - Normal BAM File
    val(norm_bai) - Normal BAI File
    val(tumor_prefix) - Tumor Prefix
    val(tumor_bam) - Tumor BAM File
    path(tumor_bai) - Tumor BAI File
  tuple
    path(fa) - Reference FASTA
    path(idx_files) - Reference Index Files
  path bed - BED File
  val parstr - Additional Parameters
  val tmp_dir - Temporary Directory

output:
  tuple => emit: norm_abra_bams
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path("*${norm_prefix}*norm_abra.bam") - Normal Realigned Output File
  tuple => emit: tumor_abra_bams
    val(pat_name) - Patient Name
    val(tumor_prefix) - Tumor Prefix
    val(norm_prefix) - Normal Prefix
    val(dataset) - Dataset
    path("*${tumor_prefix}*tumor_abra.bam") - Tumor Realigned Output File
  tuple => emit: all_abra_bams
    val(pat_name) - Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path("*${norm_prefix}*norm_abra.bam") - Normal Realigned Output File
    path("*${tumor_prefix}*tumor_abra.bam") - Tumor Realigned Output File

command:
java -jar /abra2.jar --in ${norm_bam},${tumor_bam} --out \${out_norm_bam},\${out_tumor_bam} --ref ${fa} --targets ${bed} --tmpdir ${tmp_dir} ${parstr} > abra.log
```

# abra2_cadabra

Runs abra2.jar abra.cadabra.Cadabra
Cadabra is a somatic indel caller that works specifically with ABRA alignments.
```
input:
  tuple
    val(pat_name) - Patient Name
    val(datset) - Dataset
    val(norm_prefix) - Normal Prefix
    val(norm_bam) - Normal BAM File
    val(norm_bai) - Normal BAI File
    val(tumor_prefix) - Tumor Prefix
    val(tumor_bam) - Tumor BAM File
    path(tumor_bai) - Tumor BAI File
  path fa - Reference FASTA File
  val suffix - Suffix for Output File
  val parstr - Additional Parameters

output:
  tuple => emit: vcfs
    val(pat_name) -  Patient Name
    val(norm_prefix) - Normal Prefix
    val(tumor_prefix) - Tumor Prefix
    val(dataset) - Dataset
    path("*.vcf") - Output VCF File

command:
java -cp /abra2.jar abra.cadabra.Cadabra --ref ${fa} --normal ${norm_bam} --tumor ${tumor_bam} ${parstr} > ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.${suffix}.vcf
```

# abra2_rna

Runs abra2 on RNA BAMs
ABRA is a realigner for next generation sequencing data. It uses localized
assembly and global realignment to align reads more accurately, thus
improving downstream analysis (detection of indels and complex variants in
particular).
```
input:
  tuple
    val(pat_name) - Patient Name
    val(datset) - Dataset
    val(prefix) - RNA Prefix
    val(bam) - RNA BAM File
    val(bai) - RNA BAI File
   tuple
     path(fa) - Reference FASTA File
     path(idx_files) - Reference Index Files
   path gtf - Reference GTF File
   val parstr - Additional Parameters
   val tmp_dir - Temporary Directory

 output:
   tuple => emit: abra_bams
     val(pat_name) - Patient Name
     val(norm_prefix) - Prefix
     val(dataset) - Dataset
     path("${norm_id}*abra.bam") - RNA Realigned Output File

command:
java -jar /abra2.jar --in ${bam} --out \${out_bam} --ref ${fa} ${parstr} --gtf ${gtf} --tmpdir ${tmp_dir} ${parstr} > abra.log
```

# bams_bais_to_realigned_bams

Produces ABRA2 realigned BAMs given reference FASTA, BAMs, BAIs, and manifest channels.

NOTE: There should be a 1-to-1-to-1 mapping among samples, BAMs, and BAIs
among the manifest, bams, and bais channels.

```
take:
  ref - Reference FASTA File
  bams - Channel of BAMs (see abra_rna's "abra_bams" emitted channel)
  bais - Channel of BAIs (see samtool_index's "bais" emitted channel)
  manifest - RAFT Manifest

emit:
   norm_abra_bams - Realigned BAMs from normal samples
   tumor_abra_bams - Realgined BAMS from tumor samples
   all_abra_bams - Realigned BAMs from both normal and tumor samples
```
