#hlaprofiler

HLAProfiler Nextflow Module

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| hlaprofiler_predict | docker://benjaminvincentlab/hlaprofiler:1.10.2 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |

## Workflows

`hlaprofiler_predict`

 Runs HLAProfiler.pl predict
 Uses the database included with the Docker image.
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(fq1) - FASTQ 1
    path(fq2) - FASTQ 2
  val parstr - Additional Parameters

output:
  tuple => emit: calls
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path('*') - Output Files
```
