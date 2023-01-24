# Htslib

HTSlib Nextflow Module

## Defaults

| workflow | container | cpus | memory |
| --- | --- | --- | --- |
| htslib_bgzip | docker://tedbrookings/htslib:1.9 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| htslib_bgzip_somatic | docker://tedbrookings/htslib:1.9 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |
| htslib_tabix | docker://tedbrookings/htslib:1.9 | 4 * task.attempt | 8.GB.plus(4.GB * task.attempt) |

## Workflows

`htslib_bgzip_somatic`

Runs bgzip on somatically tagged channels.
```
input:
  path inf -  Input File
  val out_dir - Output Directory

output:
  tuple => emit: tbx_file
    path(inf) -  Input File
    path("${finf}.gz") - Compressed Input File
```

`htslib_bgzip`

Runs bgzip
```
input:
  path inf -  Input File
  val out_dir - Output Directory

output:
  tuple => emit: tbx_file
    path(inf) -  Input File
    path("${inf}.gz") - Compressed Input File
```

`htslib_tabix`

Runs tabix
```
input:
  path inf -  Input File
  val out_dir - Output Directory

output:
  tuple => emit: tbx_file
    path(inf) -  Input File
    path("${finf}.gz") - Compressed Input File
```
