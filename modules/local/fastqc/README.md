# fastqc

FastQC Nextflow Module

## Defautls

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| fastqc | docker://pegi3s/fastqc:0.11.7 | 4 * task.attempt | 20.GB.plus(12.GB * task.attempt) |

## Workflows

`fastqc`

A quality control tool for high throughput sequence data.
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(fq1) - Fastq 1
    path(fq2) - Fastq 2
  val(parstr) - Additional Parameters

output:
  tuple => emit: fastqc_reports
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path("*") - Outputs
  tuple => emit: fastqc_zips
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path("*zip") - Output Zip Files
  tuple => emit: fastqc_htmls
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path("*htmls") - Output HTML Files
```
