// Load modules
params.output_dir = null
include { fastqc ; multiqc } from '../modules/Fastqc/fastqc' addParams(output_dir)

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
