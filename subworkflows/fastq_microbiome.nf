// Load modules
include { runkraken ; krakenresults ; dlkraken ; runkraken_merged ; runkraken_se} from '../modules/Microbiome/kraken2.nf' 

/*
 * Helper function to resolve Kraken database path
 * Priority: params.kraken_db > default_db_path > download
 * Returns: [db_path, needs_download]
 */
def resolveKrakenDb() {
    def default_db_path = "$baseDir/data/kraken_db/k2_minusb_20250714/"
    
    if (params.kraken_db && file(params.kraken_db).exists()) {
        log.info "[Kraken] Using user-specified database: ${params.kraken_db}"
        return [params.kraken_db, false]
    } else if (file(default_db_path).exists()) {
        log.info "[Kraken] Using default database: ${default_db_path}"
        return [default_db_path, false]
    } else {
        log.info "[Kraken] No database found - will download k2_minusb_20250714"
        return [null, true]
    }
}

/*
 * Subworkflow to get or download Kraken database
 * Emits a channel with the database path
 */
workflow KRAKEN_DB {
    main:
        def (db_path, needs_download) = resolveKrakenDb()
        
        if (needs_download) {
            dlkraken()
            kraken_db_ch = dlkraken.out.kraken_db
        } else {
            kraken_db_ch = Channel.value(file(db_path, type: 'dir'))
        }
    
    emit:
        db = kraken_db_ch
}

workflow FASTQ_KRAKEN_WF {
    take: 
        read_pairs_ch
    
    main:
        KRAKEN_DB()
        runkraken(read_pairs_ch, KRAKEN_DB.out.db)
        krakenresults(runkraken.out.kraken_report.collect())
    
    emit:
        kraken_raw = runkraken.out.kraken_raw
        kraken_report = runkraken.out.kraken_report
        krakenkrona = runkraken.out.krakenkrona_filtered
}

workflow MERGED_FASTQ_KRAKEN_WF {
    take:
        merged_reads_ch   
    
    main:
        KRAKEN_DB()
        runkraken_merged(merged_reads_ch, KRAKEN_DB.out.db)
        
        // Collect all reports (merged + unmerged)
        kraken_reports_list = runkraken_merged.out.kraken_report_merged
            .mix(runkraken_merged.out.kraken_report_unmerged)
            .collect()
        
        krakenresults(kraken_reports_list)
    
    emit:
        kraken_report_merged = runkraken_merged.out.kraken_report_merged
        kraken_report_unmerged = runkraken_merged.out.kraken_report_unmerged
}

workflow FASTQ_KRAKEN_SE_WF {
    take:
        se_nonhost_ch
    
    main:
        KRAKEN_DB()
        runkraken_se(se_nonhost_ch, KRAKEN_DB.out.db)
        krakenresults(runkraken_se.out.kraken_report.collect())
    
    emit:
        kraken_raw = runkraken_se.out.kraken_raw
        kraken_report = runkraken_se.out.kraken_report
}
