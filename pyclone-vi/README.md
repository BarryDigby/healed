# pyclone-vi

PyClone-VI Nextflow Module

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| pyclonevi_fit | docker://hub.ncsa.illinois.edu/phyloflow/pyclone-vi:latest | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| pyclonevi_write_results_file | docker://hub.ncsa.illinois.edu/phyloflow/pyclone-vi:latest | 1 * task.attempt | 2.GB.plus(2.GB * task.attempt) |

## Workflows

`pyclonevi_fit`

Runs pyclone-vi fit
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

`pyclonevi_write_results_file`

Runs pyclone-vi write_results_file
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
