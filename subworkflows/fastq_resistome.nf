// Load modules
include { index ; bwa_align ; bwa_merged_align ; samtools_merge_bams ;  samtools_merge_bams as  samtools_merge_bams_dedup} from '../modules/Alignment/bwa'

// resistome
include {plotrarefaction ; runresistome ; runsnp ; resistomeresults ; runrarefaction ; build_dependencies ; snpresults} from '../modules/Resistome/resistome'

// Deduped resistome
include { BAM_DEDUP_RESISTOME_WF } from '../subworkflows/bam_deduped_resistome.nf'

include { bwa_align_se ; samtools_dedup_se }             from '../modules/Alignment/bwa'

import java.nio.file.Paths

workflow FASTQ_RESISTOME_WF {
    take: 
        read_pairs_ch
        amr
        annotation

    main:
        // download resistome and rarefactionanalyzer
        if (file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").isEmpty()){
            build_dependencies()
            resistomeanalyzer = build_dependencies.out.resistomeanalyzer
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp =  build_dependencies.out.amrsnp
        }
        else {
            amrsnp = file("${baseDir}/bin/AmrPlusPlus_SNP/*")
            resistomeanalyzer = file("${baseDir}/bin/resistome")
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
        }
        // Define amr_index_files variable
        if (params.amr_index == null) {
            index(amr)
            amr_index_files = index.out
        } else {
            amr_index_files = Channel
                .fromPath(Paths.get(params.amr_index))
                .map { file(it.toString()) }
                .filter { file(it).exists() }
                .toList()
                .map { files ->
                    if (files.size() < 6) {
                        error "Expected 6 AMR index files, found ${files.size()}. Please provide all 6 files, including the AMR database fasta file. Remember to use * in your path."
                    } else {
                        files.sort()
                    }
                }
         }        
        // AMR alignment
        bwa_align(amr_index_files, read_pairs_ch )
        // Split sections below for standard and dedup_ed results
        runresistome(bwa_align.out.bwa_bam,amr, annotation, resistomeanalyzer )
        resistomeresults(runresistome.out.resistome_counts.collect(), "AMR")
        runrarefaction(bwa_align.out.bwa_bam, annotation, amr, rarefactionanalyzer)
        plotrarefaction(runrarefaction.out.rarefaction.collect(), "AMR")
        // Add SNP confirmation
        if (params.snp == "Y") {
            runsnp(bwa_align.out.bwa_bam, resistomeresults.out.snp_count_matrix)
            snpresults(runsnp.out.snp_counts.collect() ,"AMR" )
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
        /* ------------ (1)  DEPENDENCIES -------------------------------------- */

        if( !new File("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists() ) {
            build_dependencies()
            resistomeanalyzer   = build_dependencies.out.resistomeanalyzer
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp              = build_dependencies.out.amrsnp
        } else {
            resistomeanalyzer   = file("${baseDir}/bin/resistome")
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
            amrsnp              = file("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }
        /* ------------ (2)  AMR INDEX ----------------------------------------- */
        // declare the channel handle first

        if( params.amr_index == null ) {

            // run the indexing process you already defined
            index( amr )
            amr_index_files = index.out          // <-- channel of 6 files

        } else {

            // read files matching the user-supplied pattern
            amr_index_files = Channel
                .fromPath( params.amr_index )    // emit each path
                .ifEmpty {
                    error "No files matched '${params.amr_index}'. " +
                        "Did you forget the * wildcard?"
                }
                .collect()                       // gather into a single list
                .map { files ->
                    if( files.size() != 6 )
                        throw new RuntimeException(
                            "Expected 6 AMR index files, found ${files.size()}"
                        )
                    files.sort()                 // ensure deterministic order
                }
        }

        /* ------------ (3)  ALIGN READS --------------------------------------- */
        bwa_merged_align( amr_index_files, merged_reads_ch )

        /* ------------ (4)  MERGE BAMs ---------------------------------------- */
        def bam_pairs_ch = bwa_merged_align.out.merged_bam \
                            .mix( bwa_merged_align.out.unmerged_bam ) \
                            .groupTuple()          // (id, [bam1,bam2])

        samtools_merge_bams( bam_pairs_ch )
        def combo_bam_ch = samtools_merge_bams.out.combo_bam

        /* ------------ (5)  RESISTOME / RAREFACTION --------------------------- */
        runresistome   ( combo_bam_ch, amr, annotation, resistomeanalyzer )
        resistomeresults( runresistome.out.resistome_counts.collect() , "AMR")

        runrarefaction ( combo_bam_ch, annotation, amr, rarefactionanalyzer )
        plotrarefaction( runrarefaction.out.rarefaction.collect(), "AMR" )

        /* ------------ (6)  SNP (optional) ------------------------------------ */
        if( params.snp == 'Y' ) {
            runsnp    ( combo_bam_ch, resistomeresults.out.snp_count_matrix )
            snpresults( runsnp.out.snp_counts.collect() ,"AMR")
        }

        /* ------------ (7)  DEDUP (optional) ---------------------------------- */
        if( params.deduped == 'Y' ) {
            def dedup_pairs_ch = bwa_merged_align.out.merged_dedup_bam \
                                    .mix( bwa_merged_align.out.unmerged_dedup_bam ) \
                                    .groupTuple()
            samtools_merge_bams_dedup( dedup_pairs_ch )
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
        if( !file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists() ) {
            build_dependencies()
        }

        def amr_index_files = params.amr_index
            ? Channel.fromPath(params.amr_index)
                     .ifEmpty { error "No files matched --amr_index '${params.amr_index}'" }
                     .collect()
                     .map { it.sort() }
            : { index(amr); index.out }()

        // Align SE → coordinate-sorted BAM + index
        bwa_align_se( amr_index_files, se_nonhost_ch )

        // Optional de-dup using samtools markdup
        samtools_dedup_se( bwa_align_se.out.bwa_bam )
        def bam_for_resistome = (params.deduped == 'Y') \
            ? samtools_dedup_se.out.dedup_bam \
            : bwa_align_se.out.bwa_bam

        runresistome   ( bam_for_resistome, amr, annotation, file("${baseDir}/bin/resistome") )
        resistomeresults( runresistome.out.resistome_counts.collect(), "AMR" )

        runrarefaction ( bam_for_resistome, annotation, amr, file("${baseDir}/bin/rarefaction") )
        plotrarefaction( runrarefaction.out.rarefaction.collect(), "AMR" )

        if( params.snp == 'Y' ) {
            runsnp    ( bam_for_resistome, resistomeresults.out.snp_count_matrix )
            snpresults( runsnp.out.snp_counts.collect(), "AMR" )
        }
}