include { FASTQ_QC_WF } from "$baseDir/subworkflows/fastq_information.nf"
include { FASTQ_TRIM_WF } from './subworkflows/fastq_QC_trimming.nf'
include { FASTQ_RM_HOST_WF } from './subworkflows/fastq_host_removal.nf' 
include { FASTQ_RESISTOME_WF } from "$baseDir/subworkflows/fastq_resistome.nf"
include { FASTQ_KRAKEN_WF } from "$baseDir/subworkflows/fastq_microbiome.nf"

workflow STANDARD_AMRplusplus_wKraken {
    take: 
        read_pairs_ch
        hostfasta
        amr
        annotation
        krakendb

    main:
        // fastqc
        FASTQ_QC_WF( read_pairs_ch )
        // runqc trimming
        FASTQ_TRIM_WF(read_pairs_ch)
        // remove host DNA
        FASTQ_RM_HOST_WF(hostfasta, FASTQ_TRIM_WF.out.trimmed_reads)
        // AMR alignment
        FASTQ_RESISTOME_WF(FASTQ_RM_HOST_WF.out.nonhost_reads, amr,annotation)
        // Microbiome
        FASTQ_KRAKEN_WF(FASTQ_RM_HOST_WF.out.nonhost_reads, params.kraken_db)


}
