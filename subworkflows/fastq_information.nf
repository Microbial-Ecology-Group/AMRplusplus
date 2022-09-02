// Load modules
include { fastqc ; multiqc } from '../modules/Fastqc/fastqc' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { runqc } from '../modules/Trimming/trimmomatic' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")

// fastQC
workflow FASTQ_QC_WF {
    take: 
        read_pairs_ch

    main:
        fastqc( read_pairs_ch )
        multiqc(fastqc.out.collect(), params.multiqc )


    emit:
        fastqc = fastqc.out   
        multiqc = multiqc.out
}
