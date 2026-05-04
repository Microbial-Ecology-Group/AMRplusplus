// Switched to runnning Clumpify by default. to switch back, highlight everything under this 
// and use ctrl+F and replace to switch between "Seqkit" and "Clumpify"

// Load modules
include { SE_DeduplicateReadsClumpify} from '../modules/QC/dedup'
include { PE_DeduplicateReadsClumpify} from '../modules/QC/dedup'
include { PE_DeduplicateMergedReadsClumpify} from '../modules/QC/dedup'

// Single-end deduplication
workflow FASTQ_DEDUP_SE_WF {
  take: read_se_ch
  main:
    SE_DeduplicateReadsClumpify(read_se_ch)

  emit:
    deduped_reads = SE_DeduplicateReadsClumpify.out.dedup_fq

}

// Paired end trimming
workflow FASTQ_DEDUP_PE_WF {
  take: fastq_files
  main:
    PE_DeduplicateReadsClumpify(fastq_files)

  emit:
    deduped_reads = PE_DeduplicateReadsClumpify.out.dedup_pe_fq

}

// Paired end trimming
workflow FASTQ_DEDUP_MERGED_WF {
  take: fastq_files
  main:
    PE_DeduplicateMergedReadsClumpify(fastq_files)

  emit:
    deduped_merged_reads = PE_DeduplicateMergedReadsClumpify.out.dedup_merged_fq

}


