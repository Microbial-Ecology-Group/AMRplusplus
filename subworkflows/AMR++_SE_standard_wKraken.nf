include { FASTQ_QC_SE_WF } from "$baseDir/subworkflows/fastq_information.nf"
include { FASTQ_TRIM_SE_WF } from "$baseDir/subworkflows/fastq_QC_trimming.nf"
include { FASTQ_RM_HOST_SE_WF   } from "$baseDir/subworkflows/fastq_host_removal.nf"
include { FASTQ_RESISTOME_SE_WF } from "$baseDir/subworkflows/fastq_resistome.nf"
include { FASTQ_KRAKEN_SE_WF    } from "$baseDir/subworkflows/fastq_microbiome.nf"

workflow SE_AMRplusplus_wKraken {
    take: 
        read_se_ch
        hostfasta
        amr
        annotation

    main:
        // QC
        FASTQ_QC_SE_WF(read_se_ch)

        // Trimming
        FASTQ_TRIM_SE_WF(read_se_ch)

        // 2) Host removal (SE) on trimmed reads
        FASTQ_RM_HOST_SE_WF( hostfasta, FASTQ_TRIM_SE_WF.out.trimmed_reads )

        // 4) Branch A: Resistome
        FASTQ_RESISTOME_SE_WF( FASTQ_RM_HOST_SE_WF.out.nonhost_reads, amr, annotation )

        // 5) Branch B: Kraken microbiome
        FASTQ_KRAKEN_SE_WF( FASTQ_RM_HOST_SE_WF.out.nonhost_reads )

}


