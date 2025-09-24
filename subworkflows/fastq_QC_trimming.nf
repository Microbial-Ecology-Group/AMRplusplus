// Load modules
include { runqc ; QCstats ; QCstats_SE} from '../modules/Trimming/trimmomatic'
// Load modules (SE versions)
include { runqc_se  } from '../modules/Trimming/trimmomatic'
include { SeqkitReadCounts } from '../modules/QC/merge'

// WC trimming
workflow FASTQ_TRIM_WF {
    take: 
        read_pairs_ch

    main:
        //index( hostindex )
        //bwa_align( index.out, read_pairs_ch )
        runqc(read_pairs_ch)
        QCstats(runqc.out.trimmomatic_stats.collect())

        /* merged + unmerged QC FASTQs  → one channel  */
        runqc.out.paired_fastq
            .mix( runqc.out.unpaired_fastq )
            .set { only_qc_reads_ch }            // ← make it top-level
        
  
        /* 3 ─ one-shot SeqKit on all QC FASTQs ----------------------- */
        seqkit_input_ch = only_qc_reads_ch.map{ sid,f -> f }.collect()
        
        SeqkitReadCounts( seqkit_input_ch , "QC_trimmed" )
        
    emit:
        //bwa_align = bwa_align.out
        trimmed_reads = runqc.out.paired_fastq  
        
}


// Single-end trimming
workflow FASTQ_TRIM_SE_WF {
  take: read_se_ch
  main:
    runqc_se(read_se_ch)
    seqkit_input_ch = read_se_ch.map{ sid,f -> f }.collect()
    SeqkitReadCounts( seqkit_input_ch , "SE_QC_trimmed" )
    //QCstats_SE(runqc.out.trimmomatic_summary.collect())   // use the summary files
  emit:
    trimmed_reads = runqc_se.out.se_fastq
    //trim_stats    = QCstats_SE.out.combo_trim_stats
}