process ABRA2 {
    tag "$meta2.id"
    label 'process_high'

    conda "bioconda::abra2=2.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abra2:2.24--h7d875b9_0':
        'quay.io/biocontainers/abra2:2.24--h7d875b9_0' }"

    input:
    tuple val(meta),  path(normal_bam), path(normal_bai) // [channel] normal bam,bai or ifEmpty([])
    tuple val(meta2), path(tumor_bam),  path(tumor_bai)  // [channel] tumor  bam.bai or ifEmpty([])
    path fasta                                           // [channel] fasta
    path gtf                                             // [channel] gtf
    path junctions                                       // [channel] STAR SJ.out.tab file
    path targets                                         // [channel] intervals
    val toggle_rna                                       // [boolean] toggle RNA mode

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
    // dont use def to avail of meta information
    output_list = tumor_bam && normal_bam != [] ? "--out ${meta.id}.abra.bam,${meta2.id}.abra.bam" : tumor_bam == [] ? "--out ${meta.id}.abra.bam" : normal_bam == [] ? "--out ${meta2.id}.abra.bam" : ""
    // RNA or DNA mode means different params
    def dna_files = toggle_rna ? "" : targets ? "--targets $targets" : ""
    def rna_files = toggle_rna ? "--gtf $gtf --junctions $junctions" : ""
    def rna_args  = toggle_rna ? "--sua --dist 5000" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[ABRA2 realignment] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    java "-Xmx${avail_mem}g" -jar /usr/local/share/abra2-2.24-0/abra2.jar \\
        $input_list \\
        $output_list \\
        --tmpdir ./ \\
        --ref $fasta \\
        --threads $task.cpus \\
        $dna_files \\
        $rna_files \\
        --index \\
        $rna_args \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: \$(echo \$(abra2 2>&1 | grep 'version' | cut -d' ' -f 8))
    END_VERSIONS
    """
}
