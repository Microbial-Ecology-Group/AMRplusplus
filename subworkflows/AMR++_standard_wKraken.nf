include { FASTQ_QC_WF } from './fastq_information.nf'
include { FASTQ_TRIM_WF } from './fastq_QC_trimming.nf'
include { FASTQ_RM_HOST_WF } from './fastq_host_removal.nf' 
include { FASTQ_RESISTOME_WF } from './fastq_resistome.nf'
include { FASTQ_KRAKEN_WF } from './fastq_microbiome.nf'

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

        if( params.read_dedup == "Y" ) {
            FASTQ_DEDUP_PE_WF( FASTQ_TRIM_WF.out.trimmed_reads )
            reads_for_host = FASTQ_DEDUP_PE_WF.out.deduped_reads
        } else {
            reads_for_host = FASTQ_TRIM_WF.out.trimmed_reads
        }
        // remove host DNA
        FASTQ_RM_HOST_WF(hostfasta, reads_for_host)
        // AMR alignment
        FASTQ_RESISTOME_WF(FASTQ_RM_HOST_WF.out.nonhost_reads, amr,annotation)
        // Microbiome
        FASTQ_KRAKEN_WF(FASTQ_RM_HOST_WF.out.nonhost_reads)


}
