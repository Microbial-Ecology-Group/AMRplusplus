include { FASTQ_QC_SE_WF } from './fastq_information.nf'
include { FASTQ_TRIM_SE_WF } from './fastq_QC_trimming.nf'
include { FASTQ_RM_HOST_SE_WF   } from './fastq_host_removal.nf'
include { FASTQ_RESISTOME_SE_WF } from './fastq_resistome.nf'
include { FASTQ_KRAKEN_SE_WF    } from './fastq_microbiome.nf'
include { FASTQ_DEDUP_SE_WF } from './fastq_deduplicate.nf'

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

        if( params.read_dedup == "Y" ) {
            FASTQ_DEDUP_SE_WF( FASTQ_TRIM_SE_WF.out.se_fastq )
            reads_for_host = FASTQ_DEDUP_SE_WF.out.deduped_reads
        } else {
            reads_for_host = FASTQ_TRIM_SE_WF.out.trimmed_reads
        }

        // 2) Host removal (SE) on trimmed reads
        FASTQ_RM_HOST_SE_WF( hostfasta, reads_for_host)

        // 4) Branch A: Resistome
        FASTQ_RESISTOME_SE_WF( FASTQ_RM_HOST_SE_WF.out.nonhost_reads, amr, annotation )

        // 5) Branch B: Kraken microbiome
        FASTQ_KRAKEN_SE_WF( FASTQ_RM_HOST_SE_WF.out.nonhost_reads )

}


