include { FASTQ_QC_WF } from "$baseDir/subworkflows/fastq_information.nf"
include { FASTQ_TRIM_WF } from "$baseDir/subworkflows/fastq_QC_trimming.nf"
include { FASTQ_MERGE_WF } from "$baseDir/subworkflows/fastq_merging.nf"
include { MERGED_FASTQ_RM_HOST_WF } from "$baseDir/subworkflows/fastq_host_removal.nf" 
include { MERGED_FASTQ_RESISTOME_WF } from "$baseDir/subworkflows/fastq_resistome.nf"
include { MERGED_FASTQ_KRAKEN_WF } from "$baseDir/subworkflows/fastq_microbiome.nf"


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
        // merge reads
        FASTQ_MERGE_WF( FASTQ_TRIM_WF.out.trimmed_reads )
        
        FASTQ_MERGE_WF.out.merged
              .join( FASTQ_MERGE_WF.out.unmerged )
              .set { merged_reads_ch }
        // remove host DNA
        MERGED_FASTQ_RM_HOST_WF(hostfasta, merged_reads_ch)

        MERGED_FASTQ_RM_HOST_WF.out.nonhost_reads.into { nonhost_for_amr; nonhost_for_kraken }

        // AMR alignment
        MERGED_FASTQ_RESISTOME_WF(nonhost_for_amr, amr,annotation)

        MERGED_FASTQ_KRAKEN_WF(nonhost_for_kraken)

}


