// Load modules
include { index ; bwa_align ; bwa_merged_align ;bwa_align_se ; samtools_dedup_se ; samtools_merge_bams ;  samtools_merge_bams as  samtools_merge_bams_dedup} from '../modules/Alignment/bwa'

// resistome
include { snp_coverage_summary ; plotrarefaction ; runresistome_analyzer ; runsnp ; resistomeresults ; runrarefaction ; build_dependencies ; snpresults} from '../modules/Resistome/resistome'

// Deduped resistome
include { BAM_DEDUP_RESISTOME_WF } from '../subworkflows/bam_deduped_resistome.nf'
include { RESOLVE_BWA_INDEX } from './fastq_align.nf'


workflow FASTQ_RESISTOME_WF {
    take: 
        read_pairs_ch
        amr
        annotation

    main:
        /* ------------ (1) DEPENDENCIES ---------------------------------- */
        if ( !file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists() ) {
            build_dependencies()
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp              = build_dependencies.out.amrsnp
        } else {
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
            amrsnp              = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }

        /* ------------ (2) AMR INDEX --------------------------------------- */
        RESOLVE_BWA_INDEX(amr)
        amr_index_files = RESOLVE_BWA_INDEX.out.amr_index_files  

        /* ------------ (3) AMR ALIGNMENT ----------------------------------- */
        bwa_align(amr_index_files, read_pairs_ch )
        // Split sections below for standard and dedup_ed results
        runresistome_analyzer(bwa_align.out.bwa_bam)
        resistomeresults(runresistome_analyzer.out.resistome_counts.collect(), "AMR")
        if (params.rarefaction == "Y") {
            runrarefaction(bwa_align.out.bwa_bam, annotation, amr, rarefactionanalyzer)
            plotrarefaction(runrarefaction.out.rarefaction.collect(), "AMR")
        }
        // Add SNP confirmation
        if (params.snp == "Y") {
            runsnp(bwa_align.out.bwa_bam, resistomeresults.out.snp_count_matrix,amrsnp)
            snpresults(runsnp.out.snp_counts.collect() ,"AMR" )
            snp_coverage_summary(runsnp.out.coverage_stats.collect(), "AMR") 
        }
        // Add analysis of deduped counts
        if (params.deduped == "Y"){
            BAM_DEDUP_RESISTOME_WF(bwa_align.out.bwa_dedup_bam,amr, annotation)
        }
}


/*─────────────────────────────────────────────────────────────────────────────
 *  MERGED_FASTQ_RESISTOME_WF
 *  – work on FLASH-merged + unmerged reads
 *  – handles optional --deduped Y  and --snp Y flags
 *───────────────────────────────────────────────────────────────────────────*/
workflow MERGED_FASTQ_RESISTOME_WF {

    /* ------------ INPUTS -------------------------------------------------- */
    take:
        merged_reads_ch      // tuple(id, merged_fq, unmerged_fq)
        amr
        annotation

    main:
        /* ------------ (1) DEPENDENCIES ---------------------------------- */
        if ( !file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists() ) {
            build_dependencies()
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp              = build_dependencies.out.amrsnp
        } else {
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
            amrsnp              = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }  

        /* ------------ (2) AMR INDEX ------------------------------------- */
        RESOLVE_BWA_INDEX(amr)
        amr_index_files = RESOLVE_BWA_INDEX.out.amr_index_files  

        /* ------------ (3)  ALIGN READS --------------------------------------- */
        bwa_merged_align( amr_index_files, merged_reads_ch )

        /* ------------ (4)  MERGE BAMs ---------------------------------------- */
        def bam_pairs_ch = bwa_merged_align.out.merged_bam \
                            .mix( bwa_merged_align.out.unmerged_bam ) \
                            .groupTuple()          // (id, [bam1,bam2])

        samtools_merge_bams( bam_pairs_ch ,"")
        def combo_bam_ch = samtools_merge_bams.out.combo_bam

        /* ------------ (5)  RESISTOME / RAREFACTION --------------------------- */
        runresistome_analyzer   ( combo_bam_ch )
        resistomeresults( runresistome_analyzer.out.resistome_counts.collect() , "AMR")

        runrarefaction ( combo_bam_ch, annotation, amr, rarefactionanalyzer )
        plotrarefaction( runrarefaction.out.rarefaction.collect(), "AMR" )

        /* ------------ (6)  SNP (optional) ------------------------------------ */
        if( params.snp == 'Y' ) {
            runsnp    ( combo_bam_ch, resistomeresults.out.snp_count_matrix , amrsnp)
            snpresults( runsnp.out.snp_counts.collect() ,"AMR")
            snp_coverage_summary(runsnp.out.coverage_stats.collect(), "AMR") 
        }

        /* ------------ (7)  DEDUP (optional) ---------------------------------- */
        if( params.deduped == 'Y' ) {
            def dedup_pairs_ch = bwa_merged_align.out.merged_dedup_bam \
                                    .mix( bwa_merged_align.out.unmerged_dedup_bam ) \
                                    .groupTuple()
            samtools_merge_bams_dedup( dedup_pairs_ch , "Deduped")
            BAM_DEDUP_RESISTOME_WF( samtools_merge_bams_dedup.out.combo_bam,
                                    amr, annotation )
        }
}


workflow FASTQ_RESISTOME_SE_WF {
    take:
        se_nonhost_ch
        amr
        annotation

    main:
        /* ------------ (1) DEPENDENCIES ---------------------------------- */
        if ( !file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists() ) {
            build_dependencies()
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp              = build_dependencies.out.amrsnp
        } else {
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
            amrsnp              = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }

        /* ------------ (2) AMR INDEX --------------------------------------- */
        RESOLVE_BWA_INDEX(amr)
        amr_index_files = RESOLVE_BWA_INDEX.out.amr_index_files  

        /* ------------ (3) ALIGN SE → coordinate-sorted BAM + index --------- */
        bwa_align_se( amr_index_files, se_nonhost_ch )

        // Optional de-dup using samtools markdup
        samtools_dedup_se( bwa_align_se.out.bwa_bam )
        def bam_for_resistome = (params.deduped == 'Y') \
            ? samtools_dedup_se.out.dedup_bam \
            : bwa_align_se.out.bwa_bam

        runresistome_analyzer   ( bam_for_resistome )
        resistomeresults( runresistome_analyzer.out.resistome_counts.collect(), "AMR" )
        
        if (params.rarefaction == "Y") {
            runrarefaction ( bam_for_resistome, annotation, amr, rarefactionanalyzer )
            plotrarefaction( runrarefaction.out.rarefaction.collect(), "AMR" )
        }

        if( params.snp == 'Y' ) {
            runsnp    ( bam_for_resistome, resistomeresults.out.snp_count_matrix,amrsnp )
            snpresults( runsnp.out.snp_counts.collect(), "AMR" )
            snp_coverage_summary(runsnp.out.coverage_stats.collect(), "AMR") 
        }
}