// Load modules
include { index ; bwa_align ; bwa_merged_align } from '../modules/Alignment/bwa'

// resistome
include {plotrarefaction ; runresistome ; runsnp ; resistomeresults ; runrarefaction ; build_dependencies ; snpresults} from '../modules/Resistome/resistome'

// Deduped resistome
include { BAM_DEDUP_RESISTOME_WF } from '../subworkflows/bam_deduped_resistome.nf'

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
        resistomeresults(runresistome.out.resistome_counts.collect())
        runrarefaction(bwa_align.out.bwa_bam, annotation, amr, rarefactionanalyzer)
        plotrarefaction(runrarefaction.out.rarefaction.collect())
        // Add SNP confirmation
        if (params.snp == "Y") {
            runsnp(bwa_align.out.bwa_bam, resistomeresults.out.snp_count_matrix)
            snpresults(runsnp.out.snp_counts.collect() )
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

    take:
        merged_reads_ch    // tuple(sample_id, merged_fastq, unmerged_fastq)
        amr
        annotation

    /* ---- (1) Build or fetch auxiliary scripts -------------------- */
    Channel.create { amrsnp; resistomeanalyzer; rarefactionanalyzer }

    if( file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").isEmpty() ) {
        build_dependencies()
        amrsnp              = build_dependencies.out.amrsnp
        resistomeanalyzer   = build_dependencies.out.resistomeanalyzer
        rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
    } else {
        amrsnp              = file("${baseDir}/bin/AmrPlusPlus_SNP/*")
        resistomeanalyzer   = file("${baseDir}/bin/resistome")
        rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
    }

    /* ---- (2) Build / load AMR BWA index -------------------------- */
    if( params.amr_index == null ) {
        index(amr)
        amr_index_files = index.out
    } else {
        amr_index_files = Channel
            .fromPath(Paths.get(params.amr_index))
            .map{ file(it.toString()) }
            .filter{ it.exists() }
            .collect()
            .map { files ->
                if( files.size() < 6 ) error "Expected 6 AMR index files, found ${files.size()}"
                files.sort()
            }
    }

    /* ---- (3) Align merged + unmerged reads ----------------------- */
    bwa_merged_align( amr_index_files, merged_reads_ch )

    // Rename sample_id to keep streams separate: sampleX_merged / sampleX_unmerged
    merged_bam_ch   = bwa_merged_align.out.merged_bam   .map{ id,b -> tuple("${id}_merged",   b) }
    unmerged_bam_ch = bwa_merged_align.out.unmerged_bam .map{ id,b -> tuple("${id}_unmerged", b) }

    bam_all_ch = merged_bam_ch.mix(unmerged_bam_ch)

    /* ---- (4) Resistome + rarefaction ----------------------------- */
    runresistome   ( bam_all_ch, amr, annotation, resistomeanalyzer )
    resistomeresults( runresistome.out.resistome_counts.collect() )

    runrarefaction ( bam_all_ch, annotation, amr, rarefactionanalyzer )
    plotrarefaction( runrarefaction.out.rarefaction.collect() )

    /* ---- (5) Optional SNP verification --------------------------- */
    if( params.snp == 'Y' ) {
        runsnp     ( bam_all_ch, resistomeresults.out.snp_count_matrix )
        snpresults ( runsnp.out.snp_counts.collect() )
    }

    /* ---- (6) Optional deduped‑BAM branch ------------------------- */
    if( params.deduped == 'Y' ) {
        merged_dedup_ch   = bwa_merged_align.out.merged_dedup_bam   .map{ id,b -> tuple("${id}_merged",   b) }
        unmerged_dedup_ch = bwa_merged_align.out.unmerged_dedup_bam .map{ id,b -> tuple("${id}_unmerged", b) }

        BAM_DEDUP_RESISTOME_WF( merged_dedup_ch.mix(unmerged_dedup_ch), amr, annotation )
    }
}
