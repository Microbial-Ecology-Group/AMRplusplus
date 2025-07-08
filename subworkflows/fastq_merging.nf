// Load modules
include { index ; bwa_align } from '../modules/Alignment/bwa'

import java.nio.file.Paths

workflow FASTQ_MERGED_WF {
    take: 
        read_pairs_ch

    main:
        MergeReadsFlash( read_pairs_ch )

    emit:
        merged   = MergeReadsFlash.out.merged
        unmerged = MergeReadsFlash.out.unmerged

}


