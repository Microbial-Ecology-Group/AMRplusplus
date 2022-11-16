// Load modules
include { runkraken ; krakenresults } from '../modules/Microbiome/kraken2.nf' 

// WC trimming
workflow FASTQ_KRAKEN_WF {
    take: 
        read_pairs_ch
        krakendb

    main:
        //index( hostindex )
        //bwa_align( index.out, read_pairs_ch )
        runkraken(read_pairs_ch, krakendb)
        krakenresults(runkraken.out.kraken_report.collect())
        
}

