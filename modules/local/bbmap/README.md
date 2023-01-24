# bwa

bwa Nextflow Module

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| bwa_index | docker://fredhutch/bwa:0.7.17 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| bwa_mem | docker://fredhutch/bwa:0.7.17 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |
| bwa_mem_samtools_view | docker://fredhutch/bwa:0.7.17 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |

## Workflows

`bwa_index`

Runs bwa index
```
input:
  path fa - Reference FASTA
  val params - Additional Parameters
  val out_dir - Output Directory
  val shared_dir - Shared Output Directory

output:
  tuple => emit: idx_files
    path(fa) - Reference FASTA
    path("${fa}.*") - Index Files
```

`bwa_mem`

Runs bwa mem
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(fq1) - FASTQ 1
    path(fq1) - FASTQ 2
  tuple
    path(fa) - Reference FASTA
    path(idx_files) - Index Files
  parstr - Additional Parameters

output:
  tuple => emit: sams
    val(pat_name) - Patient Name
    val(prefix) - FASTQ prefix
    val(dataset) - Dataset
    path('*.sam') - Alignment Output File
```

`bwa_mem_samtools_view`

Runs bwa mem | samtools view -bS

NOTE: This is *not* compatible with AWS Batch.
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(fq1) - FASTQ 1
    path(fq1) - FASTQ 2
  tuple
    path(fa) - Reference FASTA
    path(idx_files) - Index Files
  parstr - Additional Parameters

output:
  tuple => emit: sams
    val(pat_name) - Patient Name
    val(prefix) - FASTQ prefix
    val(dataset) - Dataset
    path('*.sam') - Alignment Output File
```
