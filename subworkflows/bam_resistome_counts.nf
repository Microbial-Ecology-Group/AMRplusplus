// resistome counts and snp verification
include {plotrarefaction ; runresistome ; runsnp ; resistomeresults ; build_dependencies ; snpresults} from '../modules/Resistome/resistome'


workflow BAM_RESISTOME_COUNTS_WF {
    take: 
        bam_ch
        amr
        annotation

    main:
        // download resistome and rarefactionanalyzer
        if (file("${baseDir}/bin/AmrPlusPlus_SNP/SNP_Verification.py").isEmpty()){
            build_dependencies()
            resistomeanalyzer = build_dependencies.out.resistomeanalyzer
            amrsnp =  build_dependencies.out.amrsnp
        }
        else {
            amrsnp = file("${baseDir}/bin/AmrPlusPlus_SNP/*")
            resistomeanalyzer = file("${baseDir}/bin/resistome")
        }
        // Run resistome analyzer and count matrix creation
        runresistome(bam_ch,amr, annotation, resistomeanalyzer )
        resistomeresults(runresistome.out.resistome_counts.collect())
        // Add SNP confirmation
        if (params.snp == "Y") {
            runsnp(bam_ch, resistomeresults.out.snp_count_matrix)
            snpresults(runsnp.out.snp_counts.collect(), "AMR")
        }
}


