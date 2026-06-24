// Load modules
include { index ; bwa_align ; bwa_merged_align ;bwa_align_se ; samtools_dedup_se ; samtools_merge_bams ;  samtools_merge_bams as  samtools_merge_bams_dedup } from '../modules/Alignment/bwa'
include { build_dependencies } from '../modules/Resistome/resistome'

workflow FASTQ_ALIGN_WF {
    take: 
        read_pairs_ch
        amr

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
            amrsnp              = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
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


workflow MERGED_FASTQ_ALIGN_WF {

    /* ------------ INPUTS -------------------------------------------------- */
    take:
        merged_reads_ch      // tuple(id, merged_fq, unmerged_fq)
        amr

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
            amrsnp              = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }

        /* ------------ (2) AMR INDEX ------------------------------------- */
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

        /* ------------ (3)  ALIGN READS --------------------------------------- */
        bwa_merged_align( amr_index_files, merged_reads_ch )

        /* ------------ (4)  MERGE BAMs ---------------------------------------- */
        def bam_pairs_ch = bwa_merged_align.out.merged_bam \
                            .mix( bwa_merged_align.out.unmerged_bam ) \
                            .groupTuple()          // (id, [bam1,bam2])

        samtools_merge_bams( bam_pairs_ch )
}

workflow SE_FASTQ_ALIGN_WF {
    take:
        se_nonhost_ch
        amr

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
            amrsnp              = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
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

        /* ------------ (3) ALIGN SE → coordinate-sorted BAM + index --------- */
        bwa_align_se( amr_index_files, se_nonhost_ch )

        // Optional de-dup using samtools markdup
        samtools_dedup_se( bwa_align_se.out.bwa_bam )
        def bam_for_resistome = (params.deduped == 'Y') \
            ? samtools_dedup_se.out.dedup_bam \
            : bwa_align_se.out.bwa_bam
}