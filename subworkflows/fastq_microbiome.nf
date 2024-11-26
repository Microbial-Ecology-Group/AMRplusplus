// Load modules
include { runkraken ; krakenresults ; dlkraken} from '../modules/Microbiome/kraken2.nf' 

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
