// Load modules
include { runkraken ; krakenresults ; dlkraken} from '../modules/Microbiome/kraken2.nf' 

// WC trimming
workflow FASTQ_KRAKEN_WF {
    take: 
        read_pairs_ch
        krakendb

    main:
        if (params.kraken_db == null) {
            dlkraken()
            kraken_db_ch = dlkraken.out
        } else {
            kraken_db_ch = Channel
                .fromPath(params.kraken_db)
         }    
        runkraken(read_pairs_ch, kraken_db_ch)
        krakenresults(runkraken.out.kraken_report.collect())
        
}

