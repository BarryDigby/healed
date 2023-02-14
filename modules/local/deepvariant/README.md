# deepvariant

DeepVariant Nextflow Module

## Defautls

| workflow    | container                         | cpus              | memory                            |
| ----------- | --------------------------------- | ----------------- | --------------------------------- |
| deepvariant | docker://google/deepvariant:1.1.0 | 4 \* task.attempt | 20.GB.plus(12.GB \* task.attempt) |

## Workflows

`deepvariant`

Runs DeepVariant on a sample

```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    path(bam) - Alignment File
  path(fa) - Reference FASTA
  model_type - Model Type
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
