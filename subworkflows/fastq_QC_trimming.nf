// Load modules
include { runqc ; QCstats ; QCstats_SE} from '../modules/Trimming/trimmomatic'
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
  take: read_se_ch
  main:
    runqc_se(read_se_ch)
    //QCstats_SE(runqc.out.trimmomatic_summary.collect())   // use the summary files
  emit:
    trimmed_reads = runqc_se.out.se_fastq
    //trim_stats    = QCstats_SE.out.combo_trim_stats
}