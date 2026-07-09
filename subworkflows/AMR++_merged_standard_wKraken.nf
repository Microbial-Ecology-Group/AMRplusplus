include { FASTQ_QC_WF } from './fastq_information.nf'
include { FASTQ_TRIM_WF } from './fastq_QC_trimming.nf'
include { FASTQ_MERGE_WF } from './fastq_merging.nf'
include { MERGED_FASTQ_RM_HOST_WF } from './fastq_host_removal.nf' 
include { MERGED_FASTQ_RESISTOME_WF } from './fastq_resistome.nf'
include { MERGED_FASTQ_KRAKEN_WF } from './fastq_microbiome.nf'
include { FASTQ_DEDUP_PE_WF } from './fastq_deduplicate.nf'

workflow STANDARD_merged_AMRplusplus_wKraken {
    take: 
        read_pairs_ch
        hostfasta
        amr
        annotation

    main:
        // fastqc
        FASTQ_QC_WF( read_pairs_ch )
        // runqc trimming
        FASTQ_TRIM_WF(read_pairs_ch)
        
        if( params.read_dedup == "Y" ) {
            FASTQ_DEDUP_PE_WF( FASTQ_TRIM_WF.out.trimmed_reads )
            reads_for_merge = FASTQ_DEDUP_PE_WF.out.deduped_reads
        } else {
            reads_for_merge = FASTQ_TRIM_WF.out.trimmed_reads
        }
        // merge reads
        FASTQ_MERGE_WF( reads_for_merge )
        
        FASTQ_MERGE_WF.out.merged
              .join( FASTQ_MERGE_WF.out.unmerged )
              .set { merged_reads_ch }
        // remove host DNA
        MERGED_FASTQ_RM_HOST_WF(hostfasta, merged_reads_ch)

        // AMR alignment
        MERGED_FASTQ_RESISTOME_WF(MERGED_FASTQ_RM_HOST_WF.out.nonhost_reads, amr,annotation)

        // Microbiome kraken2
        MERGED_FASTQ_KRAKEN_WF(MERGED_FASTQ_RM_HOST_WF.out.nonhost_reads)

}


