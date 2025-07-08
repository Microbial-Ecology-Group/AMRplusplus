// Load modules
include { index } from '../modules/Alignment/bwa'
include { bwa_align ; bwa_rm_contaminant_fq ; bwa_rm_contaminant_merged_fq; HostRemovalStats} from '../modules/Alignment/bwa'
include { SeqkitReadCounts } from '../modules/QC/merge'



import java.nio.file.Paths

//  Host removal for paired reads
workflow FASTQ_RM_HOST_WF {
    take: 
        hostfasta
        read_pairs_ch
    main:
        // Define reference_index variable
        if (params.host_index == null) {
            index(hostfasta)
            reference_index_files = index.out
        } else {
            reference_index_files = Channel
                .fromPath(Paths.get(params.host_index))
                .map { file(it.toString()) }
                .filter { file(it).exists() }
                .toList()
                .map { files ->
                    if (files.size() < 6) {
                        error "Expected 6 host index files, found ${files.size()}. Please provide all 6 files, including the host fasta file. Remember to use * in your path."
                    } else {
                        files.sort()
                    }
                }
         }    
        bwa_rm_contaminant_fq(reference_index_files, read_pairs_ch )
        HostRemovalStats(bwa_rm_contaminant_fq.out.host_rm_stats.collect())
    emit:
        nonhost_reads = bwa_rm_contaminant_fq.out.nonhost_reads  
}

workflow MERGED_FASTQ_RM_HOST_WF {
    take: 
        hostfasta
        merged_reads_ch

    main:
        /* 1 ─ build / load BWA index -------------------------------------- */
        def reference_index_ch =
            params.host_index
            ? Channel.fromPath( params.host_index , glob:true )
                    .ifEmpty { error "No files match --host_index '${params.host_index}'" }
                    .toList()
                    .map { it.sort() }               // bundle 6 index files
            : { index( hostfasta ); index.out }()     // call in a closure


        /* 2 ─ host-removal -------------------------------------------------- */
        bwa_rm_contaminant_merged_fq( reference_index_ch , merged_reads_ch )
        
        bwa_rm_contaminant_merged_fq.out.nonhost_merged
            .join( bwa_rm_contaminant_merged_fq.out.nonhost_unmerged )
            .set { nonhost_reads_ch } 
        
        /* merged + unmerged non-host FASTQs  → one channel  */
        bwa_rm_contaminant_merged_fq.out.nonhost_merged
            .mix( bwa_rm_contaminant_merged_fq.out.nonhost_unmerged )
            .set { only_reads_ch }            // ← make it top-level
        
        
        /* 3 ─ one-shot SeqKit on all non-host FASTQs ----------------------- */
        seqkit_input_ch = only_reads_ch.map{ sid,f -> f }.collect()
        
        SeqkitReadCounts( seqkit_input_ch , "NonHost" )



    emit:
        nonhost_reads_ch
}