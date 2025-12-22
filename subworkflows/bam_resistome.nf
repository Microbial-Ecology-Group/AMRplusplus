// resistome
include {temp_runsnp; plotrarefaction ; runresistome ; runsnp ; resistomeresults ; runrarefaction ; build_dependencies ; snpresults} from '../modules/Resistome/resistome'


workflow BAM_RESISTOME_WF {
    take: 
        bam_ch
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
        // Split sections below for standard and dedup_ed results
        runresistome(bam_ch,amr, annotation, resistomeanalyzer )
        resistomeresults(runresistome.out.resistome_counts.collect(), "AMR")
        // Rarefaction
        if (params.rarefaction == "Y") {
            runrarefaction(bam_ch, annotation, amr, rarefactionanalyzer)
            plotrarefaction(runrarefaction.out.rarefaction.collect(), "AMR")
        }
        // Add SNP confirmation
        if (params.snp == "Y") {
            temp_runsnp(bam_ch, resistomeresults.out.snp_count_matrix,file("${baseDir}/bin/AmrPlusPlus_SNP/*"))
            snpresults(runsnp.out.snp_counts.collect(), "AMR")
        }
}


