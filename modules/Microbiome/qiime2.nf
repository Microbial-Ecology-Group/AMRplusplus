p_trim_left_f = params.p_trim_left_f
p_trim_left_r = params.p_trim_left_r
p_trunc_len_f = params.p_trunc_len_f
p_trunc_len_r = params.p_trunc_len_r

threads = params.threads

process Qiime2Import {
    tag { }
    label "small"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Qiime2Results", mode: "copy"

    input:
        path(manifest)

    output:
        path("demux.qza"), emit: demux

    """
    ${QIIME} tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path ${manifest} \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33       
    """
}


process Qiime2Dada2 {
    tag { }
    label "medium"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Qiime2Results", mode: "copy"

    input:
        path(demux)

    output:
        path("dada-table.qza"), emit: dada_table
        path("rep-seqs.qza"), emit: rep_seqs

    """
    ${QIIME} dada2 denoise-paired --i-demultiplexed-seqs ${demux} --o-table dada-table.qza --o-representative-sequences rep-seqs.qza --p-trim-left-f ${p_trim_left_f} --p-trim-left-r ${p_trim_left_r} --p-trunc-len-f ${p_trunc_len_f} --p-trunc-len-r ${p_trunc_len_r} --p-n-threads ${threads} --verbose --o-denoising-stats denoise_stats 

    """
}



process Qiime2Classify {
    tag { }
    label "medium"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Qiime2Results", mode: "copy"

    input:
        path(rep_seqs)
        path(database)

    output:
        path("taxonomy.qza"), emit: taxonomy

    """
    ${QIIME} feature-classifier classify-sklearn --i-classifier ${database} --i-reads ${rep_seqs} --o-classification taxonomy.qza

    """
}

process Qiime2Filter {
    tag { }
    label "medium"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Qiime2Results", mode: "copy"

    input:
        path(dada_table)
        path(taxonomy)
        path(rep_seqs)

    output:
        path("filtered_table.qza"), emit: filtered_table
        path("filtered_rep-seqs.qza"), emit: filtered_seqs

    """
    ${QIIME} taxa filter-table --i-table ${dada_table} --i-taxonomy ${taxonomy} --p-exclude mitochondria,chloroplast --o-filtered-table filtered_table.qza 

    ${QIIME} taxa filter-seqs --i-sequences ${rep_seqs} --i-taxonomy ${taxonomy} --p-exclude mitochondria,chloroplast --o-filtered-sequences filtered_rep-seqs.qza

    """
}

process Qiime2Tree {
    tag { }
    label "medium"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Qiime2Results/Tree", mode: "copy"

    input:
        path(filtered_seqs)

    output:
        path("rooted-tree.qza"), emit: rooted_tree
        
    """
    ${QIIME} alignment mafft --i-sequences ${filtered_seqs} --o-alignment aligned-rep-seqs.qza 
    
    ${QIIME} alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza 
    
    ${QIIME} phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza 
    
    ${QIIME} phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza 
    """
}

process Qiime2Export {
    tag { }
    label "medium"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/Qiime2Results/Exported", mode: "copy"

    input:
        path(filtered_seqs)
        path(rooted_tree)
        path(filtered_table)
        path(taxonomy)

    output:
        path("table-with-taxonomy.biom"), emit: table_taxa
        path("tree.nwk"), emit: tree_nwk
        path("dna-sequences.fasta"), emit: dna_seqs

    """
    ${QIIME} tools export --input-path filtered_rep-seqs.qza --output-path .
    ${QIIME} tools export --input-path taxonomy.qza --output-path . 
    ${QIIME} tools export --input-path rooted-tree.qza --output-path .
    ${QIIME} tools export --input-path filtered_table.qza --output-path . 

    # Change out column headers in taxonomy file
    sed -i 's/Feature ID/#OTUID/g' taxonomy.tsv
    sed -i 's/Taxon/taxonomy/g' taxonomy.tsv
    sed -i 's/Confidence/confidence/g' taxonomy.tsv

    # Add taxonomy file to biom
    biom add-metadata -i feature-table.biom -o table-with-taxonomy.biom --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy

    rm feature-table.biom
    """
}