process ARRIBA_REF_DOWNLOAD {
    tag 'ARRIBA_download'

    conda "bioconda::arriba=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.3.0--haa8aa89_0' :
        'quay.io/biocontainers/arriba:2.3.0--haa8aa89_0' }"

    output:
    path "blacklist_hg38_GRCh38_v2.3.0.tsv.gz"            , emit: blacklist
    path "known_fusions_hg38_GRCh38_v2.3.0.tsv.gz"        , emit: known_fusions
    path "protein_domains_hg38_GRCh38_v2.3.0.gff3"        , emit: protein_domains


    script:
    """
    wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz --no-check-certificate

    tar xvf arriba_v2.3.0.tar.gz

    mv arriba_v2.3.0/database/blacklist_hg38_GRCh38_v2.3.0.tsv.gz .
    mv arriba_v2.3.0/database/protein_domains_hg38_GRCh38_v2.3.0.gff3 .
    mv arriba_v2.3.0/database/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz .

    rm -rf arriba_v2.3.0*
    """
}
