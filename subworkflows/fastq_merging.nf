// Load modules
include {MergeReadsFlash } from '../modules/QC/merge'

import java.nio.file.Paths

workflow FASTQ_MERGE_WF {
    take: 
        read_pairs_ch

    main:
        MergeReadsFlash( read_pairs_ch )

    emit:
        merged   = MergeReadsFlash.out.merged
        unmerged = MergeReadsFlash.out.unmerged

}


