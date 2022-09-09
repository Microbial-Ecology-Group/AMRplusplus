// Load modules
include { index ; bwa_align } from '../modules/Alignment/bwa'

// resistome
include { runresistome ; runsnp ; resistomeresults ; runrarefaction ; build_dependencies} from '../modules/Resistome/resistome'


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
            amrsnp = build_dependencies.out.amrsnp
            // Index
            index(amr)
            // AMR alignment
            bwa_align(amr, index.out, read_pairs_ch )
            runresistome(bwa_align.out.bwa_sam,amr, annotation, resistomeanalyzer )
            runsnp(bwa_align.out.bwa_sam, amrsnp )
            resistomeresults(runresistome.out.resistome_counts.collect())
            runrarefaction(bwa_align.out.bwa_sam, annotation, amr, rarefactionanalyzer)
        }
        else {
            amrsnp = file("${baseDir}/bin/AmrPlusPlus_SNP/")
            resistomeanalyzer = file("${baseDir}/bin/resistome")
            rarefactionanalyzer = file("${baseDir}/bin/rarefaction")
            // Index
            index(amr)
            // AMR alignment
            bwa_align(amr, index.out, read_pairs_ch )
            runresistome(bwa_align.out.bwa_sam,amr, annotation, resistomeanalyzer )
            runsnp(bwa_align.out.bwa_sam, amrsnp )
            resistomeresults(runresistome.out.resistome_counts.collect())
            runrarefaction(bwa_align.out.bwa_sam, annotation, amr, rarefactionanalyzer)
        }


}
