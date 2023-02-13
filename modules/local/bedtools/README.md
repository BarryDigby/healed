# Bedtools

BEDTools Nextflow Module

Website: https://bedtools.readthedocs.io/en/latest/

## Defaults

| workflow           | container                                   | cpus              | memory                          |
| ------------------ | ------------------------------------------- | ----------------- | ------------------------------- |
| bedtools_intersect | docker://biocontainers/bedtools:v2.28.0_cv2 | 4 \* task.attempt | 8.GB.plus(8.GB \* task.attempt) |

## Workflows

# bedtools_intersect

Compares two or more BED/BAM/VCF/GFF files and identifies all the regions in
the genome where the features in the two files overlap (that is, share at
least one base pair in common).

```
input:
  tuple
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path(aln) - Alignment File
 path bed - BED file

output:
  tuple => emit: filt_alns
    val(pat_name) - Patient Name
    val(prefix) - Prefix
    val(dataset) - Dataset
    path("*filt.bam") - BAM file reduced to BED file coordinates

command:
bedtools intersect -a ${aln} -b ${bed}  > \${NEW_NAME}.filt.bam
```
