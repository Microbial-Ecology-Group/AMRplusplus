// resistome counts and snp verification
include {plotrarefaction ; runresistome_analyzer ; runsnp ; resistomeresults ; build_dependencies ; snpresults} from '../modules/Resistome/resistome'


workflow BAM_RESISTOME_COUNTS_WF {
    take: 
        bam_ch
        amr
        annotation

    main:
        // download resistome and rarefactionanalyzer
        if (file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").isEmpty()){
            build_dependencies()
            amrsnp =  build_dependencies.out.amrsnp
        }
        else {
            amrsnp = file("${baseDir}/bin/AmrPlusPlus_SNP/*")
        }
        // Run resistome analyzer and count matrix creation
        runresistome_analyzer(bam_ch )
        resistomeresults(runresistome.out.resistome_counts.collect(), "AMR")
        // Add SNP confirmation
        if (params.snp == "Y") {
            runsnp(bam_ch, resistomeresults.out.snp_count_matrix)
            snpresults(runsnp.out.snp_counts.collect(), "AMR")
        }
}


