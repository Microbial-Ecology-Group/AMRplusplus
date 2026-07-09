include { FASTQ_QC_WF} from './fastq_information.nf'
include { FASTQ_TRIM_WF } from './fastq_QC_trimming.nf'
include { FASTQ_RESISTOME_WF } from './fastq_resistome.nf'

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
        
        if( params.read_dedup == "Y" ) {
            FASTQ_DEDUP_PE_WF( FASTQ_TRIM_WF.out.trimmed_reads )
            reads_for_host = FASTQ_DEDUP_PE_WF.out.deduped_reads
        } else {
            reads_for_host = FASTQ_TRIM_WF.out.trimmed_reads
        }
        // AMR alignment
        FASTQ_RESISTOME_WF(reads_for_host, amr,annotation)

        
}
