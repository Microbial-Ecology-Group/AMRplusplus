// Load modules
include { runqc ; QCstats } from '../modules/Trimming/trimmomatic'
// Load modules (SE versions)
include { runqc_se  } from '../modules/Trimming/trimmomatic'

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


// Single-end trimming
workflow FASTQ_TRIM_SE_WF {
    take:
        read_se_ch   // emits: tuple(sample_id, read.fastq[.gz])

    main:
        runqc_se(read_se_ch)
        QCstats(runqc_se.out.trimmomatic_stats.collect())

    emit:
        trimmed_reads = runqc_se.out.se_fastq
        trim_stats    = QCstats.out.combo_trim_stats
}