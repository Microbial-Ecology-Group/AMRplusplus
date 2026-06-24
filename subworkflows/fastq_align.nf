// Load modules
include { index ; bwa_align } from '../modules/Alignment/bwa'

import java.nio.file.Paths

workflow FASTQ_ALIGN_WF {
    take: 
        read_pairs_ch
        amr
        annotation

    main:
        /* ------------ (1) DEPENDENCIES ---------------------------------- */
        if ( !file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists() ) {
            build_dependencies()
            resistomeanalyzer   = build_dependencies.out.resistomeanalyzer
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp              = build_dependencies.out.amrsnp
        } else {
            resistomeanalyzer   = file("${baseDir}/bin/resistome")
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
            amrsnp              = file("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }

        /* ------------ (2) AMR INDEX --------------------------------------- */
        if (params.amr_index) {
            amr_index_files = Channel
                .fromPath(params.amr_index, glob: true)
                .ifEmpty { error "No files match --amr_index '${params.amr_index}'" }
                .collect()
                .map { files ->
                    if (files.size() < 7) {
                        error "Expected 7 AMR index files, found ${files.size()}. Please provide all 7 files, including the AMR database fasta file. Remember to use * in your path."
                    }
                    files.sort()
                }
        } else {
            index(amr)
            amr_index_files = index.out
        }     
        /* ------------ (3) AMR ALIGNMENT ----------------------------------- */
        bwa_align(amr_index_files, read_pairs_ch )
}


