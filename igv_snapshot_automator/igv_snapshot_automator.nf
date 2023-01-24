process igv_snapshot_automator {

    tag "${dataset}/${pat_name}"

    label 'igv_snapshot_automator_container'
    label 'igv_snapshot_automator'

    publishDir "${params.lens_out_dir}/${dataset}/${pat_name}/igv_snapshot_automator"

    input:
    tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai), val(rna_run), path(rna_bam), path(rna_bai), path(bed)

    output:
    tuple val(pat_name), val(dataset), path("*igv"), path("*igv.bams")

    script:
    """
    python /IGV-snapshot-automator/make_IGV_snapshots.py \
      -r ${bed} \
      -g hg38 \
      -o \${PWD}/${dataset}-${pat_name}-igv \
      -nf4 \
      ${norm_bam} \
      ${tumor_bam} \
      ${rna_bam} \
      -ht 5000

    mkdir -p ${dataset}-${pat_name}.igv.bams
    mv *bam ${dataset}-${pat_name}.igv.bams
    """
}
