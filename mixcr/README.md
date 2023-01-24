# mixcr

MiXCR Nextflow Module

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| mixcr_align | docker://benjaminvincentlab/mixcr:3.0.13.1 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| mixcr_assemble | docker://benjaminvincentlab/mixcr:3.0.13.1 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| mixcr_export | docker://benjaminvincentlab/mixcr:3.0.13.1 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| mixcr_shotgun | docker://benjaminvincentlab/mixcr:3.0.13.1 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |

## Workflows

`mixcr_align`

Runs mixcr align
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
  tuple => emit: vdjca_files
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path('*vdjca') - VDJCA Output Files
```

`mixcr_assemble`

Runs mixcr assemble
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(vdjca) - VDJCA Input File
  val if_config - null = assemble, non-null = assembleContigs
  val parstr - Additional Parameters

output:
  tuple => emit: cln_files
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path('*.cln*') - Output File (either *.clna or *.clns)
```

`mixcr_export`

Runs mixcr export
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path(infile) - Input File (either *.clna or *.clns)
  val exports - "clones" = exportClones, *.clna and "alignments" = exportAlignments
  val parstr - Additional Parameters

output:
  tuple => emit: exported_files
    val(pat_name) - Patient Name
    val(prefix) - FASTQ Prefix
    val(dataset) - Dataset
    path('*txt') - Output file
```

`mixcr_shotgun`

Runs mixcr analyze shotgun & exportReads
```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(datasest) - Dataset
    path(fq1) - FASTQ 1
    path(fq2) - FASTQ2
  val(regex) - Output Regular Expression
  val(parstr) - Additional Parameters

output:
  tuple => emit: clonotypes
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(datasest) - Dataset
    path("${regex}") - exportClones Clonotypes
  tuple => emit: report
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(datasest) - Dataset
    path("*.report") - Combined Report
  tuple => emit: aligned_r1
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(datasest) - Dataset
    path("*_aligned_r1.fq.gz" - mixcr align R1 Reads
  tuple => emit: aligned_r2
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(datasest) - Dataset
    path("*_aligned_r2.fq.gz" - mixcr align R2 Reads
```
