// Load modules
include { index as amr_index ; index as host_index } from '../modules/Alignment/bwa'
include { bwa_align ; bwa_rm_contaminant_fq ; HostRemovalStats} from '../modules/Alignment/bwa'


if( params.host_index ) {
    host_index_files = Channel.fromPath(params.reference_index).toSortedList()
    //if( host_index.isEmpty() ) return index_error(host_index)
}

// WC trimming
workflow FASTQ_RM_HOST_WF {
    take: 
        hostfasta
        read_pairs_ch
    main:
        if( !params.host_index_files ) {    
            host_index(hostfasta)
            host_index_files = host_index.out.bwaindex
        }
        bwa_rm_contaminant_fq(hostfasta,host_index_files, read_pairs_ch )
        HostRemovalStats(bwa_rm_contaminant_fq.out.host_rm_stats.collect())
    emit:
        nonhost_reads = bwa_rm_contaminant_fq.out.nonhost_reads  
}

