include { FASTQ_QC_WF } from "$baseDir/subworkflows/fastq_information.nf"
include { FASTQ_TRIM_WF } from "$baseDir/subworkflows/fastq_QC_trimming.nf"
include { FASTQ_MERGE_WF } from "$baseDir/subworkflows/fastq_merging.nf"
include { MERGED_FASTQ_RM_HOST_WF } from "$baseDir/subworkflows/fastq_host_removal.nf" 
include { MERGED_FASTQ_RESISTOME_WF } from "$baseDir/subworkflows/fastq_resistome.nf"

workflow STANDARD_merged_AMRplusplus {
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

        // AMR alignment
        MERGED_FASTQ_RESISTOME_WF(MERGED_FASTQ_RM_HOST_WF.out.nonhost_reads, amr,annotation)


}


