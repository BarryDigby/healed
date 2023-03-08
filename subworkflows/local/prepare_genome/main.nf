//
// PREPARE GENOME
//

// DNA
include { BWA_INDEX                                       } from '../../../modules/nf-core/bwa/index/main'

// RNA-Seq
include { STAR_GENOMEGENERATE                             } from '../../../modules/nf-core/star/genomegenerate/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/rsem/preparereference/main'
include { GTF_GENE_FILTER                                 } from '../../../modules/local/gtf_gene_filter/main'

//RNA Fusion
include { ARRIBA_REF_DOWNLOAD                             } from '../../../modules/local/arriba/ref_download/main'
include { CTAT_GENOME_LIB_DOWNLOAD                        } from '../../../modules/local/ctat_genome/download/main'
inlcude { CTAT_GENOME_LIB_BUILD                           } from '../../../modules/local/ctat_genome/build/main'
include { FUSIONCATCHER_DOWNLOAD                          } from '../../../modules/local/fusioncatcher/download/main'

workflow PREPARE_GENOME {
    take:
        fasta                   // channel: [mandatory] fasta
        fasta_fai               // channel: [optional]  fasta_fai
        gtf                     // channel: [mandatory] gtf
        prepare_tool_indices    // [Array.list]


    main:

    ch_versions = Channel.empty()

    ch_arriba_blacklist       = Channel.empty()
    ch_arriba_cytobands       = Channel.empty()
    ch_arriba_known_fusions   = Channel.empty()
    ch_arriba_protein_domains = Channel.empty()
    if('arriba' in prepare_tool_indices) {
        if( (params.arriba_blacklist && params.arriba_cytobands && params.arriba_known_fusions && params.arriba_protein_domains) ) {
            ch_arriba_blacklist       = file(params.arriba_blacklist)
            ch_arriba_cytobands       = file(params.arriba_cytobands)
            ch_arriba_known_fusions   = file(params.arriba_known_fusions)
            ch_arriba_protein_domains = file(params.arriba_protein_domains)
        } else {
        ARRIBA_REF_DOWNLOAD()
        ch_arriba_blacklist       = ARRIBA_REF_DOWNLOAD.out.blacklist
        ch_arriba_cytobands       = ARRIBA_REF_DOWNLOAD.out.cytobands
        ch_arriba_known_fusions   = ARRIBA_REF_DOWNLOAD.out.known_fusions
        ch_arriba_protein_domains = ARRIBA_REF_DOWNLOAD.out.protein_domains
        }
    }

    ch_bwa_index = Channel.empty()
    if('bwa' in prepare_tool_indices) {
        if(params.bwa_index) {
            ch_bwa_index = Channel.fromPath(params.bwa_index).collect()
            ch_bwa_index = ch_bwa_index.map{ it -> [[id:it[0].baseName], it]}
        } else {
            ch_bwa_index = BWA_INDEX(fasta.map{ it -> [[id:it[0].baseName], it] }).index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    ch_star_index = Channel.empty()
    //if('star' in prepare_tool_indices)
    if(prepare_tool_indices.every{["star", "star_fusion", "arriba"]}) {
        if (params.star_index) {
            ch_star_index = Channel.fromPath(params.star_index).collect()
        } else {
            ch_star_index = STAR_GENOMEGENERATE( fasta, gtf ).index
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    ch_ctat_genome_lib = Channel.empty()
    ch_ctat_genome_gtf = Channel.empty()
    if('star_fusion' in prepare_tool_indices) {
        if (params.ctat_genome_lib) {
            ch_ctat_genome_lib = file(params.ctat_genome_lib)
            ch_ctat_genome_gtf = file("${params.ctat_genome_lib}/ref_annot.gtf")
        } else {
            CTAT_GENOME_LIB_DOWNLOAD()
            ch_ctat_genome_lib = CTAT_GENOME_LIB_DOWNLOAD.out.reference
            ch_ctat_genome_gtf = CTAT_GENOME_LIB_DOWNLOAD.out.chrgtf
        }
    }

    ch_filter_gtf       = Channel.empty()
    ch_transcript_fasta = Channel.empty()
    if('star' in prepare_tool_indices) {
        if(params.transcript_fasta){
            ch_transcript_fasta = file(params.transcript_fasta)
        } else {
            ch_filter_gtf       = GTF_GENE_FILTER ( fasta, gtf ).gtf
            ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( fasta, ch_filter_gtf ).transcript_fasta
            ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
            ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
        }
    }

    ch_fusioncatcher_lib = Channel.empty()
    if('fusioncatcher' in prepare_tool_indices) {
        if(params.fusioncatcher_lib){
            ch_fusioncatcher_lib = Channel.fromPath(params.fusioncatcher_lib)
        } else {
            FUSIONCATCHER_DOWNLOAD()
            ch_fusioncatcher_lib = FUSIONCATCHER_DOWNLOAD.out.reference
            ch_versions          = ch_versions.mix(FUSIONCATCHER_DOWNLOAD.out.versions)
        }
    }

    // Gather versions of all tools used
    //ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    // ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions)
/*     ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_versions = ch_versions.mix(MSISENSORPRO_SCAN.out.versions)
    ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
    ch_versions = ch_versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_SNPS.out.versions)
    ch_versions = ch_versions.mix(TABIX_KNOWN_INDELS.out.versions)
    ch_versions = ch_versions.mix(TABIX_PON.out.versions) */
  //  ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
        arriba_blacklist                 = ch_arriba_blacklist
        arriba_cytobands                 = ch_arriba_cytobands
        arriba_known_fusions             = ch_arriba_known_fusions
        arriba_protein_domains           = ch_arriba_protein_domains
        bwa_index                        = ch_bwa_index       // path: bwa/*
        ctat_genome_lib                  = ch_ctat_genome_lib
        ctat_genome_gtf                  = ch_ctat_genome_gtf
/*         bwamem2                          = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwamem2/*
        hashtable                        = DRAGMAP_HASHTABLE.out.hashmap.map{ meta, index -> [index] }.collect() // path: dragmap/*
        dbsnp_tbi                        = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
        germline_resource_tbi            = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()   // path: germline_resource.vcf.gz.tbi
        known_snps_tbi                   = TABIX_KNOWN_SNPS.out.tbi.map{ meta, tbi -> [tbi] }.collect()          // path: {known_indels*}.vcf.gz.tbi
        known_indels_tbi                 = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        msisensorpro_scan                = MSISENSORPRO_SCAN.out.list.map{ meta, list -> [list] }                // path: genome_msi.list
        pon_tbi                          = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()                 // path: pon.vcf.gz.tbi
        chr_files                        = chr_files
        allele_files                     = allele_files
        loci_files                       = loci_files
        gc_file                          = gc_file
        rt_file                          = rt_file */
        filter_gtf                       = ch_filter_gtf
        fusioncatcher_lib                = ch_fusioncatcher_lib
        star_index                       = ch_star_index
        transcript_fasta                 = ch_transcript_fasta
        versions                         = ch_versions                                                         // channel: [ versions.yml ]
}
