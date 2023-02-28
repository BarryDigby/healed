process ABRA2 {
    tag "$meta2.id"
    label 'process_high'

    conda "bioconda::abra2=2.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abra2:2.24--h9f5acd7_1':
        'quay.io/biocontainers/abra2:2.24--h9f5acd7_1' }"

    input:
    tuple val(meta),  path(normal_bam), path(normal_bai) // [channel] normal bam,bai or ifEmpty([])
    tuple val(meta2), path(tumor_bam),  path(tumor_bai)  // [channel] tumor  bam.bai or ifEmpty([])
    path fasta                                           // [channel] fasta
    path gtf                                             // [channel] gtf
    path junctions                                       // [channel] STAR SJ.out.tab file
    path targets                                         // [channel] intervals

    output:
    tuple val(meta),  path("${meta.id}.abra.bam"),  path("${meta.id}.abra.bai"),    optional:true, emit: abra2_normal
    tuple val(meta2), path("${meta2.id}.abra.bam"), path("${meta2.id}.abra.bai"),   optional:true, emit: abra2_tumor
    path "versions.yml",                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // define input lists. Use .ifEmpty([]) on inputs to invoke below logic.
    def input_list = tumor_bam && normal_bam != [] ? "--in $normal_bam,$tumor_bam" : tumor_bam == [] ? "--in $normal_bam" : normal_bam == [] ? "--in $tumor_bam" : ""
    output_list = tumor_bam && normal_bam != [] ? "--out ${meta.id}.abra.bam,${meta2.id}.abra.bam" : tumor_bam == [] ? "--out ${meta.id}.abra.bam" : normal_bam == [] ? "--out ${meta2.id}.abra.bam" : ""
    def targs = targets ? "--targets $targets" : ""
    def gtf_file = gtf ? "--gtf $gtf" : ""
    def juncs = junctions ? "--junctions $junctions" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[ABRA2 realignment] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    // Not sure what conda path actually is but handler in place below
    def binPath = ( params.enable_conda ? "abra2" : "/usr/local/share/abra2-2.24-1/abra2.jar" )
    """
    java "-Xmx${avail_mem}g" -jar $binPath \\
        $input_list \\
        $output_list \\
        --tmpdir ./ \\
        --ref $fasta \\
        --threads $task.cpus \\
        $targs \\
        $gtf_file \\
        --index \\
        $juncs \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: \$(echo \$(abra2 2>&1 | grep 'version' | cut -d' ' -f 8))
    END_VERSIONS
    """
}
