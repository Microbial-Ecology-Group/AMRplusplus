// Load modules
include { index } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { bwa_align } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")

// resistome
include { runresistome } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")
include { runsnp } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")
include { resistomeresults } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")
include { runrarefaction } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")

workflow FASTQ_RESISTOME_WF {
    take: 
        read_pairs_ch
        amr
        annotation

    main:
        index(amr)
        // AMR alignment
        bwa_align(amr, index.out, read_pairs_ch )
        runresistome(bwa_align.out.bwa_sam,amr, annotation )
        //runsnp(bwa_align.out.bwa_sam )
        resistomeresults(runresistome.out.resistome_counts.collect())
        runrarefaction(bwa_align.out.bwa_sam, annotation, amr)
    emit:
        rarefaction_results = runrarefaction.out.rarefaction


}
