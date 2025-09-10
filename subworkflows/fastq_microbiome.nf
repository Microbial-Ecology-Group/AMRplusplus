// Load modules
include { runkraken ; krakenresults ; dlkraken ; runkraken_merged ; runkraken_se} from '../modules/Microbiome/kraken2.nf' 

workflow FASTQ_KRAKEN_WF {
    take: 
        read_pairs_ch
        krakendb

    main:
        // Define the default database path
        def default_db_path = "$baseDir/data/kraken_db/minikraken_8GB_20200312/"
        
        // Prioritize params.kraken_db over the default_db_path if defined
        def db_path = params.kraken_db ?: (file(default_db_path).exists() ? default_db_path : null)
        
        // Logic to handle database availability
        if (db_path == null) {
            // If no database is available, download one
            dlkraken()
            runkraken(read_pairs_ch, dlkraken.out)
        } else {
            // Use the available database (params.kraken_db or default_db_path)
            kraken_db_ch = Channel.value(db_path)
            runkraken(read_pairs_ch, kraken_db_ch)
        }
        krakenresults(runkraken.out.kraken_report.collect())
}


workflow MERGED_FASTQ_KRAKEN_WF {

    take:
        merged_reads_ch   

    main:

        
        /*  Build ( sid , mergedFile , unmergedFile )  -------------------------- */

          
          // now each item is exactly: ( sid, Path-to-merged, Path-to-unmerged )
                
        /*──────────────── choose / download Kraken DB ───────────────*/
        def kraken_db_ch
        def default_db = "$baseDir/data/kraken_db/k2_standard_08gb_20250402/"
        if( file(default_db).exists() )
             kraken_db_ch = Channel.value(default_db)
        else if( params.kraken_db )
             kraken_db_ch = Channel.value(params.kraken_db)
        else {
             dlkraken()
             kraken_db_ch = dlkraken.out
        }

        /*──────────────── Kraken 2  ─────────────────────────────────*/
        runkraken_merged( merged_reads_ch , kraken_db_ch )

        
        /* ---------- run krakenresults once all reports are ready ---------- */
        def kraken_reports_list = runkraken_merged.out.kraken_report_merged
                                    .mix( runkraken_merged.out.kraken_report_unmerged )
                                    .collect()                       // Java List (one per sample)
        
        krakenresults( kraken_reports_list )   // ← plain value, no Channel.value()
              
          
        /* ---------- one-shot SeqKit on all extracted FASTQs --------------- */
        //seqkit_fastq_list = only_extracted_reads_ch.map{ sid,f -> f }.collect()
        //SeqkitReadCounts( seqkit_fastq_list , "Kraken_extracted" )


            
            
            
}

workflow FASTQ_KRAKEN_SE_WF {
    take:
        se_nonhost_ch
        krakendb

    main:
        def default_db = "$baseDir/data/kraken_db/k2_standard_08gb_20250402/"
        def db_ch =
            (params.kraken_db ? Channel.value(params.kraken_db)
          : file(default_db).exists() ? Channel.value(default_db)
          : { dlkraken(); dlkraken.out }() )

        runkraken_se( se_nonhost_ch, db_ch )
        krakenresults( runkraken_se.out.kraken_report.collect() )
}