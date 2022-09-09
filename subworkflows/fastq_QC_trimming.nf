// Load modules
include { runqc } from '../modules/Trimming/trimmomatic'

// WC trimming
workflow FASTQ_TRIM_WF {
    take: 
        read_pairs_ch

    main:
        //index( hostindex )
        //bwa_align( index.out, read_pairs_ch )
        runqc(read_pairs_ch)
    emit:
        //bwa_align = bwa_align.out
        trimmed_reads = runqc.out.paired_fastq  
        
}

