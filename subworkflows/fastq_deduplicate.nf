// Load modules
include { DeduplicateReadsSeqkit} from '../modules/QC/dedup'



// Single-end trimming
workflow FASTQ_DEDUP_SE_WF {
  take: read_se_ch
  main:
    DeduplicateReadsSeqkit(read_se_ch)

  emit:
    deduped_reads = DeduplicateReadsSeqkit.out.dedup_fq

}



