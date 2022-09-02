// Load modules
include { index as amr_index ; index as host_index } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { bwa_align ; bwa_rm_contaminant_fq} from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")

// WC trimming
workflow FASTQ_RM_HOST_WF {
    take: 
        hostfasta
        read_pairs_ch
    main:
        host_index(hostfasta)
        bwa_rm_contaminant_fq(hostfasta, host_index.out, read_pairs_ch )
    emit:
        nonhost_reads = bwa_rm_contaminant_fq.out.nonhost_reads  
}

