// Load modules
include { index as amr_index ; index as host_index } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { bwa_align } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { fastqc ; multiqc } from '../modules/Fastqc/fastqc' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { runqc } from '../modules/Trimming/trimmomatic' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")

// WC trimming
workflow FASTQ_TRIM_WF {
    take: 
        read_pairs_ch

    main:
        //index( hostindex )
        //bwa_align( index.out, read_pairs_ch )
        
        runqc(read_pairs_ch)
        
        //multiqc(bwa_align.out.mix(fastqc.out).collect(), params.multiqc )
        //multiqc(fastqc.out.collect(), params.multiqc )

    emit:
        //bwa_align = bwa_align.out
        trimmed_reads = runqc.out.paired_fastq  
        
}
