// Load modules
include { runqc ; QCstats } from '../modules/Trimming/trimmomatic'

// WC trimming
workflow FASTQ_TRIM_WF {
    take: 
        read_pairs_ch

    main:
        //index( hostindex )
        //bwa_align( index.out, read_pairs_ch )
        runqc(read_pairs_ch)
        QCstats(runqc.out.trimmomatic_stats.collect())
        
    emit:
        //bwa_align = bwa_align.out
        trimmed_reads = runqc.out.paired_fastq  
        
}

