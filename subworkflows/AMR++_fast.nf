include { FASTQ_QC_WF } from "$baseDir/subworkflows/fastq_information.nf"
include { FASTQ_RESISTOME_WF } from "$baseDir/subworkflows/fastq_resistome.nf"

workflow FAST_AMRplusplus {
    take: 
        read_pairs_ch
        amr
        annotation

    main:
        // fastqc
        FASTQ_QC_WF( read_pairs_ch )
        // runqc trimming
        FASTQ_TRIM_WF(read_pairs_ch)
        // AMR alignment
        FASTQ_RESISTOME_WF(FASTQ_TRIM_WF.out.trimmed_reads, amr,annotation)

        
}
