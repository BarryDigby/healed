# picard2

Picard2 Nextflow Module

## Defaults

| `| container | cpus | memory |
| --- | --- | --- | --- |
| picard_create_seq_dict | docker://broadinstitute/picard:2.21.4 | 2 * task.attempt | 4.GB.plus(4.GB * task.attempt) |
| picard_mark_duplicates | docker://broadinstitute/picard:2.21.4 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| picard_collect_insert_size_metrics | docker://broadinstitute/picard:2.21.4 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| picard_collect_rna_seq_metrics | docker://broadinstitute/picard:2.21.4 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| picard_collect_wgs_metrics_nzc | docker://broadinstitute/picard:2.21.4 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| picard_collect_vcf_metrics | docker://broadinstitute/picard:2.21.4 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |

## Workflows

`picard_create_seq_dict`

Runs CreateSequenceDictionary
```
input:
  path(fa) - Reference FASTA
  val(parstr) - Additional Parameters

output:
  tuple => emit: dict_file
    path(fa) - Reference FASTA
    path('*.dict') - Dict File
```

`picard_mark_duplicates`

Runs MarkDuplicates
```
input:
  path aln - Alignment Files
  val parstr - Additional Parameters

output:
  tuple => emit: mkdup_alns
    path(fa) - Reference FASTA
    path('*.dict') - Dict File
```

`picard_collect_insert_size_metrics`

Runs CollectInsertSizeMetrics
```
input:
  path aln - Alignment Files
  val parstr - Additional Parameters

output:
  tuple => emit: insert_metrics
    path(fa) - Reference FASTA
    path('*.dict') - Dict File
```

`picard_collect_rna_seq_metrics`

Runs CollectRnaSeqMetrics
```
input:
output:
```

`picard_collect_wgs_metrics_nzc`

Collect metrics about coverage and performance of whole genome sequencing
(WGS) experiments. This tool collects metrics about the percentages of reads
that pass base- and mapping- quality filters as well as coverage
(read-depth) levels. Both minimum base- and mapping-quality values as well
as the maximum read depths (coverage cap) are user defined. This extends
CollectWgsMetrics by including metrics related only to siteswith non-zero
(>0) coverage.
```
input:
output:
```

`picard_collect_vcf_metrics`

Collects per-sample and aggregate (spanning all samples) metrics from the
provided VCF file.
```
input:
output:
```
