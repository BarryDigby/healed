# salmon

Salmon Nextflow Module

## Defaults

| `                | container                        | cpus              | memory                          |
| ---------------- | -------------------------------- | ----------------- | ------------------------------- |
| salmon_index     | docker://combinelab/salmon:1.1.0 | 4 \* task.attempt | 8.GB.plus(4.GB \* task.attempt) |
| salmon_map_quant | docker://combinelab/salmon:1.1.0 | 4 \* task.attempt | 8.GB.plus(4.GB \* task.attempt) |
| salmon_aln_quant | docker://combinelab/salmon:1.1.0 | 4 \* task.attempt | 8.GB.plus(4.GB \* task.attempt) |

## Workflows

`salmon_index`

Runs salmon index command

```
input:
  path fa - Reference FASTA
  val parstr - Additional Parameters

output:
  tuple => emit: idx_files
    path(fa) - Reference FASTA
    path("${fa}*index") - Index Files
```

`salmon_map_quant`

Runs salmon quant (for quantifying from FASTQs)

```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(fq1) - FASTQ 1
    path(fq2) - FASTQ 2
  tuple
    path(fa) - Reference FA
    path(idx_files) - Index Files
  val parstr - Additoinal Parameters

output:
  tuple => emit: quants
      val(pat_name) - Patient Name
      val(prefix) - FASTQ Prefix
      val(dataset) - Dataset
      path('*quant.sf') - Quant File
```

`salmon_aln_quant`

Runs salmon quant (for quantifying from alignments)

```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(aln) - Alignment File
  path(fa) - Reference FA
  val out_dir - Output Directory
  val shared_dir - Shared Output Directory

output:
  tuple => emit: quants
      val(pat_name) - Patient Name
      val(prefix) - FASTQ Prefix
      val(dataset) - Dataset
      path('*quant.sf') - Quant File
```
