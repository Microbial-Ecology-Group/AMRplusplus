// resistome
include {snp_coverage_summary;coverage_threshold_sweep ; combine_sweep_results ; plot_sweep_dropoff ; plotrarefaction ; runresistome_analyzer ; runsnp ; resistomeresults ; runrarefaction ; build_dependencies ; snpresults  } from '../modules/Resistome/resistome'


workflow BAM_RESISTOME_WF {
    take: 
        bam_ch
        amr
        annotation

    main:
        // download resistome and rarefactionanalyzer
        if (!file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").exists()){
            build_dependencies()
            rarefactionanalyzer = build_dependencies.out.rarefactionanalyzer
            amrsnp =  build_dependencies.out.amrsnp
        }
        else {
            amrsnp = files("${baseDir}/bin/AmrPlusPlus_SNP/*")
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
        }
        // Split sections below for rarefaction and SNP confirmation
        runresistome_analyzer(bam_ch)
        resistomeresults(runresistome_analyzer.out.resistome_counts.collect(), "AMR")
        // Rarefaction
        if (params.rarefaction == "Y") {
            runrarefaction(bam_ch, annotation, amr, rarefactionanalyzer)
            plotrarefaction(runrarefaction.out.rarefaction.collect(), "AMR")
        }
        // Add SNP confirmation
        if (params.snp == "Y") {
            runsnp(bam_ch, resistomeresults.out.snp_count_matrix)
            snpresults(runsnp.out.snp_counts.collect(), "AMR")
            snp_coverage_summary(runsnp.out.coverage_stats.collect(), "AMR") 
        }
}

workflow BAM_COVERAGE_SWEEP_WF {
    take:
        bam_ch     // tuple(prefix, bam)  -- prefix = BAM basename minus .bam

    main:
        // 1. One sweep per BAM. Each emits its four <prefix>_*.csv files.
        coverage_threshold_sweep( bam_ch )

        // 2. Gather EVERY per-sample CSV across ALL samples into one list, then
        //    combine into the four combined_*.csv matrices. flatten() because each
        //    sweep task emits a list of 4 files; collect() to wait for all samples.
        all_sweep_csvs = coverage_threshold_sweep.out.sweep_csvs.flatten().collect()
        combine_sweep_results( all_sweep_csvs )

        // 3. Plot dropoff summaries from the combined matrices.
        plot_sweep_dropoff( combine_sweep_results.out.combined )

    emit:
        combined = combine_sweep_results.out.combined
}
