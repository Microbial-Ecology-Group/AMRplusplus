// Load modules
include { runkraken ; krakenresults ; dlkraken} from '../modules/Microbiome/kraken2.nf' 

workflow FASTQ_KRAKEN_WF {
    take: 
        read_pairs_ch
        krakendb

    main:
        if (params.kraken_db == null) {
            if (file("$baseDir/data/kraken_db/minikraken_8GB_20200312/").isDirectory()) {
                kraken_db_ch = Channel.fromPath("$baseDir/data/kraken_db/minikraken_8GB_20200312/")
            } else {
                dlkraken()
                kraken_db_ch = dlkraken.out
            }
        } else {
            kraken_db_ch = Channel.fromPath(params.kraken_db)
        }
        
        runkraken(read_pairs_ch, kraken_db_ch)
        krakenresults(runkraken.out.kraken_report.collect())
        
}

