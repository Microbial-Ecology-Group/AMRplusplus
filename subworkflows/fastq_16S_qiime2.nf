
// resistome
include { Qiime2Import  ; Qiime2Dada2 ; Qiime2Classify ; Qiime2Filter ; Qiime2Tree ; Qiime2Export } from '../modules/Microbiome/qiime2'


workflow FASTQ_QIIME2_WF {
    take: 
        manifest
        database
        
    main:
        Qiime2Import(manifest)
        Qiime2Dada2(Qiime2Import.out.demux)
        Qiime2Classify(Qiime2Dada2.out.rep_seqs, database)
        Qiime2Filter(Qiime2Dada2.out.dada_table, Qiime2Classify.out.taxonomy , Qiime2Dada2.out.rep_seqs)
        Qiime2Tree(Qiime2Filter.out.filtered_seqs)
        Qiime2Export(Qiime2Filter.out.filtered_seqs, Qiime2Tree.out.rooted_tree, Qiime2Filter.out.filtered_table, Qiime2Classify.out.taxonomy )

}
