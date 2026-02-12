// Load modules
include { fastqc ; multiqc } from '../modules/Fastqc/fastqc'
include { SeqkitReadCounts } from '../modules/QC/merge'

// fastQC
workflow FASTQ_QC_WF {
    take: 
        read_pairs_ch

    main:
        fastqc( read_pairs_ch )
        multiqc(fastqc.out.collect(), params.multiqc )
        //seqkit_input_ch = read_pairs_ch.map{ sid,f -> f }.collect()
        //SeqkitReadCounts( seqkit_input_ch , "Raw" )
}

// Single-end FastQC + MultiQC
workflow FASTQ_QC_SE_WF {
    take:
        read_se_ch   // emits: tuple(sample_id, read.fastq[.gz])

    main:
        fastqc(read_se_ch)
        // pass FastQC outputs + your MultiQC config (params.multiqc) to MultiQC
        multiqc(fastqc.out.collect(), params.multiqc)

    emit:
        fastqc = fastqc.out
}