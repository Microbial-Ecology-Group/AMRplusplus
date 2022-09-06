params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
params.readlen = 150

process runkraken {
    tag { sample_id }
    conda = "$baseDir/envs/microbiome.yaml"
    container = 'enriquedoster/amrplusplus_microbiome:latest'

    publishDir "${params.output}/RunKraken", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".kraken.raw") > 0) "Standard/$filename"
            else if(filename.indexOf(".kraken.report") > 0) "Standard_report/$filename"
            else if(filename.indexOf(".kraken.filtered.report") > 0) "Filtered_report/$filename"
            else if(filename.indexOf(".kraken.filtered.raw") > 0) "Filtered/$filename"
            else {}
        }

    input:
       tuple val(sample_id), path(reads)
       path(krakendb)


   output:
      tuple val(sample_id), path("${sample_id}.kraken.raw"), emit: kraken_raw
      path("${sample_id}.kraken.report"), emit: kraken_report
      tuple val(sample_id), path("${sample_id}.kraken.filtered.raw"), emit: kraken_filter_raw
      path("${sample_id}.kraken.filtered.report"), emit: kraken_filter_report
      tuple val(sample_id), path("${sample_id}_kraken2.krona"), emit: krakenkrona
      tuple val(sample_id), path("${sample_id}_kraken2_filtered.krona"), emit: krakenkrona_filtered



     """
     ${KRAKEN2} --db ${krakendb} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
     ${KRAKEN2} --db ${krakendb} --confidence 1 --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.filtered.raw

    cut -f 2,3  ${sample_id}.kraken.raw > ${sample_id}_kraken2.krona
    cut -f 2,3  ${sample_id}.kraken.filtered.raw > ${sample_id}_kraken2_filtered.krona
    """
}

process krakenresults {
    tag { }
    conda = "$baseDir/envs/microbiome.yaml"
    container = 'enriquedoster/amrplusplus_microbiome:latest'
    publishDir "${params.output}/KrakenResults", mode: "copy"

    input:
        path(kraken_reports)

    output:
        path("kraken_analytic_matrix.csv")

    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o kraken_analytic_matrix.csv
    """
}

process runbracken {

    input:
       tuple val(sample_id), path(krakenout)
       tuple val(sample_id), path(krakenout_filtered)
       path(krakendb)

    """
    bracken \
        -d ${krakendb} \
        -r ${params.readlen} \
        -i ${krakenout} \
        -l ${params.taxlevel} \
        -o ${sample_id}_bracken.tsv

    bracken \
        -d ${krakendb} \
        -r ${params.readlen}\
        -i ${krakenout_filtered} \
        -l ${params.taxlevel} \
        -o ${sample_id}_bracken_filtered.tsv
        """
}

process kronadb {
    output:
        file("krona_db/taxonomy.tab") optional true into krona_db_ch // is this a value ch?

    when: 
        !params.skip_krona
        
    script:
    """
    ktUpdateTaxonomy.sh krona_db
    """
}

process kronafromkraken {
    publishDir params.outdir, mode: 'copy'

    input:
        file(x) from kraken2krona_ch.collect()
        //file(y) from kaiju2krona_ch.collect()
        file("krona_db/taxonomy.tab") from krona_db_ch
    
    output:
        file("*_taxonomy_krona.html")

    when:
        !params.skip_krona
    
    script:
    """
    mkdir krona
    ktImportTaxonomy -o kraken2_taxonomy_krona.html -tax krona_db $x
    """
}