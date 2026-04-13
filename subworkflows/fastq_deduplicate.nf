// Load modules
include { SE_DeduplicateReadsSeqkit} from '../modules/QC/dedup'
include { PE_DeduplicateReadsSeqkit} from '../modules/QC/dedup'
include { PE_DeduplicateMergedReadsSeqkit} from '../modules/QC/dedup'



// Single-end deduplication
workflow FASTQ_DEDUP_SE_WF {
  take: read_se_ch
  main:
    SE_DeduplicateReadsSeqkit(read_se_ch)

  emit:
    deduped_reads = SE_DeduplicateReadsSeqkit.out.dedup_fq

}

// Paired end trimming
workflow FASTQ_DEDUP_PE_WF {
  take: fastq_files
  main:
    PE_DeduplicateReadsSeqkit(fastq_files)

  emit:
    deduped_reads = PE_DeduplicateReadsSeqkit.out.dedup_pe_fq

}

// Paired end trimming
workflow FASTQ_DEDUP_MERGED_WF {
  take: fastq_files
  main:
    PE_DeduplicateMergedReadsSeqkit(fastq_files)

  emit:
    deduped_merged_reads = PE_DeduplicateMergedReadsSeqkit.out.dedup_merged_fq

}


