// Load modules
include { index ; bwa_align ; bwa_merged_align ;bwa_align_se ; samtools_dedup_se ; samtools_merge_bams ;  samtools_merge_bams as  samtools_merge_bams_dedup } from '../modules/Alignment/bwa'
include { build_dependencies } from '../modules/Resistome/resistome'

workflow RESOLVE_BWA_INDEX {

    take:
        amr          // path to the reference fasta (--amr)

    main:
        use_supplied_index = false

        if (params.amr_index) {
            amr_name  = file(amr).name
            idx_files = files(params.amr_index)
            idx_bases = idx_files.collect {
                it.name.replaceAll(/\.(amb|ann|bwt|pac|sa|fai)$/, '')
            } as Set

            if (idx_bases.isEmpty()) {
                log.warn """
                --amr_index '${params.amr_index}' matched no files. Building a new
                index from --amr instead.
                """.stripIndent()
            } else if (idx_bases.contains(amr_name)) {
                use_supplied_index = true
            } else {
                log.warn """
                ────────────────────────────────────────────────────────────────
                --amr_index does not match --amr, so it is being IGNORED and a new
                index will be built from --amr.

                  --amr           : ${amr}
                  --amr_index     : ${params.amr_index}
                  index built for : ${idx_bases.join(', ')}

                If this is not what you intended, check --amr for a typo. To use a
                pre-built index, point --amr_index at the index for ${amr_name}.
                ────────────────────────────────────────────────────────────────
                """.stripIndent()
            }
        }

        if (use_supplied_index) {
            log.info "Using pre-built BWA index: ${params.amr_index}"
            amr_index_files = Channel
                .fromPath(params.amr_index, glob: true)
                .ifEmpty { error "No files match --amr_index '${params.amr_index}'" }
                .collect()
                .map { fs ->
                    if (fs.size() < 7) {
                        error "Expected 7 AMR index files, found ${fs.size()}. " +
                              "Please provide all 7 files, including the AMR database " +
                              "fasta file. Remember to use * in your path."
                    }
                    fs.sort()
                }
        } else {
            log.info "Building BWA index from: ${amr}"
            index(amr)
            amr_index_files = index.out
        }

    emit:
        amr_index_files
}


workflow FASTQ_ALIGN_WF {
    take: 
        read_pairs_ch
        amr

    main:
        /* ------------ (1) AMR INDEX --------------------------------------- */
        RESOLVE_BWA_INDEX(amr)
        amr_index_files = RESOLVE_BWA_INDEX.out.amr_index_files  
        /* ------------ (2) AMR ALIGNMENT ----------------------------------- */
        bwa_align(amr_index_files, read_pairs_ch )
}


workflow MERGED_FASTQ_ALIGN_WF {

    /* ------------ INPUTS -------------------------------------------------- */
    take:
        merged_reads_ch      // tuple(id, merged_fq, unmerged_fq)
        amr

    main:
        /* ------------ (1) AMR INDEX ------------------------------------- */
        RESOLVE_BWA_INDEX(amr)
        amr_index_files = RESOLVE_BWA_INDEX.out.amr_index_files  

        /* ------------ (2)  ALIGN READS --------------------------------------- */
        bwa_merged_align( amr_index_files, merged_reads_ch )

        /* ------------ (3)  MERGE BAMs ---------------------------------------- */
        def bam_pairs_ch = bwa_merged_align.out.merged_bam \
                            .mix( bwa_merged_align.out.unmerged_bam ) \
                            .groupTuple()          // (id, [bam1,bam2])

        samtools_merge_bams( bam_pairs_ch,"standard" )

        // Add analysis of deduped counts
        if (params.deduped == "Y"){
            def bam_deduped_pairs_ch = bwa_merged_align.out.merged_dedup_bam \
                            .mix( bwa_merged_align.out.unmerged_dedup_bam ) \
                            .groupTuple()          // (id, [bam1,bam2])
            samtools_merge_bams_dedup ( bam_pairs_ch,"deduped" )


        }

}

workflow SE_FASTQ_ALIGN_WF {
    take:
        se_nonhost_ch
        amr

    main:
        /* ------------ (1) AMR INDEX --------------------------------------- */
        RESOLVE_BWA_INDEX(amr)
        amr_index_files = RESOLVE_BWA_INDEX.out.amr_index_files  

        /* ------------ (2) ALIGN SE → coordinate-sorted BAM + index --------- */
        bwa_align_se( amr_index_files, se_nonhost_ch )

        // Optional de-dup using samtools markdup
        samtools_dedup_se( bwa_align_se.out.bwa_bam )
        def bam_for_resistome = (params.deduped == 'Y') \
            ? samtools_dedup_se.out.dedup_bam \
            : bwa_align_se.out.bwa_bam
}