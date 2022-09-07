// Load modules
include { fastqc ; multiqc } from '../modules/Fastqc/fastqc'

// fastQC
workflow FASTQ_QC_WF {
    take: 
        read_pairs_ch

    main:
        fastqc( read_pairs_ch )
        multiqc(fastqc.out.collect(), params.multiqc )

}
