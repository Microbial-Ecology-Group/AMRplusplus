// subworkflows/bam_snv.nf
//
// SNV calling starting from BAM files, analogous in shape to BAM_RESISTOME_WF.
//
// Pipeline shape:
//   bam_ch
//     -> ngless_filter_snv   per sample: min length, min identity, unique only
//     -> index_snv_bams      per sample: metaSNV needs an index
//     -> run_metasnv         ONCE across all samples
//     -> clean_metasnv       raw calls -> AMR++-style count matrix + annotations
//     -> filter_snv_matrix   optional: restrict to genes seen in the resistome
//
// See modules/SNV/metasnv.nf for why the reference passed here must be the same
// fasta the BAMs were aligned to, and why that should normally be the MEGARes
// representative database.

include { ngless_filter_snv ; index_snv_bams ; run_metasnv ;
          clean_metasnv ; filter_snv_matrix } from '../modules/Resistome/metasnv.nf'

workflow BAM_SNV_WF {

    take:
        bam_ch             // tuple(sample_id, bam)
        reference          // fasta the BAMs were aligned to
        resistome_matrix   // AMR_analytic_matrix.csv, or an empty channel

    main:
        /* ── (1) NGLess filtering ───────────────────────────────────────────
         * Minimum aligned-block size, minimum percent identity, unique mapping
         * only. The unique-mapping requirement is the important part: metaSNV
         * cannot attribute an allele to a gene when the read maps equally well
         * elsewhere.
         */
        ngless_filter_snv( bam_ch )

        /* ── (2) index ──────────────────────────────────────────────────────
         * NGLess samtools_sort already produces coordinate-sorted output, so
         * only the index is missing.
         */
        index_snv_bams( ngless_filter_snv.out.filtered_bam )

        /* ── (3) metaSNV, once across all samples ───────────────────────────*/
        all_bams = index_snv_bams.out.indexed_bam.map { id, bam, bai -> bam }.collect()
        all_bais = index_snv_bams.out.indexed_bam.map { id, bam, bai -> bai }.collect()

        run_metasnv( all_bams, all_bais, reference )

        /* ── (4) parse into an AMR++-style matrix ───────────────────────────*/
        clean_metasnv( run_metasnv.out.snv_dir )

        /* ── (5) optional: restrict to genes present in the resistome ───────*/
        if (params.snv_filter_by_resistome == "Y") {
            filter_snv_matrix( clean_metasnv.out.snv_matrix, resistome_matrix )
        }

    emit:
        snv_matrix      = clean_metasnv.out.snv_matrix
        snv_annotations = clean_metasnv.out.snv_annotations
        snv_dir         = run_metasnv.out.snv_dir
}
