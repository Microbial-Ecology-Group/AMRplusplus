process runkraken {
    tag { sample_id }
    conda = "$baseDir/envs/AMR++_microbiome.yaml"
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


     """
     ${KRAKEN2} --db ${krakendb} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
     ${KRAKEN2} --db ${krakendb} --confidence 1 --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.filtered.raw
     """
}

process krakenresults {
    tag { }
    conda = "$baseDir/envs/AMR++_microbiome.yaml"
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