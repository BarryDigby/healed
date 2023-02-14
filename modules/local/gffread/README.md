# gffread

gffread Nextflow Module

## Defaults

| workflow           | container                       | cpus              | memory                          |
| ------------------ | ------------------------------- | ----------------- | ------------------------------- |
| gffread_make_tx_fa | docker://zavolab/gffread:0.11.7 | 1 \* task.attempt | 1.GB.plus(3.GB \* task.attempt) |

## Workflows

`gffread_make_tx_fa`
GFF/GTF utility providing format conversions, filtering, FASTA sequence
extraction and more.

```
input:
  path(fa) - Reference FASTA
  path(gtf) - Reference GTF

output:
  path("transcripts.fa") => emit: tx_fa
```
