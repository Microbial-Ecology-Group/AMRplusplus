params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
params.readlen = 150

threads = params.threads
kraken_confidence = params.kraken_confidence

process dlkraken {
    tag { }
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "$baseDir/data/kraken_db/", mode: 'copy'

    output:
        path("minikraken_8GB_20200312/")

    """
        wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
        tar -xvzf minikraken_8GB_202003.tgz

    """
}

process runkraken {
    tag { sample_id }
    label "microbiome"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/MicrobiomeAnalysis", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".conf_${kraken_confidence}.kraken.raw") > 0) "Kraken/Raw_output_conf_${kraken_confidence}/$filename"
            else if(filename.indexOf(".conf_${kraken_confidence}.kraken.report") > 0) "Kraken/Report_conf_${kraken_confidence}/$filename"
            else {}
        }

    input:
       tuple val(sample_id), path(reads)
       path(krakendb)


   output:
      tuple val(sample_id), path("${sample_id}.conf_${kraken_confidence}.kraken.raw"), emit: kraken_raw
      path("${sample_id}.conf_${kraken_confidence}.kraken.report"), emit: kraken_report
      tuple val(sample_id), path("${sample_id}.conf_${kraken_confidence}.kraken.krona"), emit: krakenkrona_filtered


     """
     ${KRAKEN2} --db ${krakendb} --confidence ${kraken_confidence} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.conf_${kraken_confidence}.kraken.report > ${sample_id}.conf_${kraken_confidence}.kraken.raw

     cut -f 2,3  ${sample_id}.conf_${kraken_confidence}.kraken.raw > ${sample_id}.conf_${kraken_confidence}.kraken.krona
    """
}

process runkraken_merged {

    tag   { sample_id }
    // robust to null | String | List
    // make sure to test this functionality
    label {
        def optsStr =
            (params.kraken_options instanceof List) ? params.kraken_options.join(' ')
          : (params.kraken_options != null)         ? params.kraken_options.toString()
                                                    : ''

        optsStr.contains('--memory-mapping') ? 'small_memory_long_time' : 'large_long'
    }

    publishDir "${params.output}/MicrobiomeAnalysis", mode: 'copy',
        saveAs: { fn ->
            if      (fn.endsWith('.kraken.raw'))   "Kraken/standard/$fn"
            else if (fn.endsWith('.kraken.report'))"Kraken/standard_report/$fn"
            else if (fn.endsWith('.fastq.gz'))     "Kraken/extracted_reads/$fn"
        }

    input:
        tuple val(sample_id), path(merged), path(unmerged)   // now BOTH are single files
        val krakendb

    output:
        tuple val(sample_id), path("${sample_id}.merged.kraken.raw"),      emit: kraken_raw_merged
        path("${sample_id}.merged.kraken.report"),                         emit: kraken_report_merged
        tuple val(sample_id), path("${sample_id}_extracted_merged.fastq.gz"), emit: extracted_merged

        tuple val(sample_id), path("${sample_id}.unmerged.kraken.raw"),    emit: kraken_raw_unmerged
        path("${sample_id}.unmerged.kraken.report"),                       emit: kraken_report_unmerged
        tuple val(sample_id), path("${sample_id}_extracted_unmerged.fastq.gz"), emit: extracted_unmerged

    script:
    """
    # ── merged file ─────────────────────────────────────────────
    ${KRAKEN2} --db ${krakendb} ${kraken_options} --confidence ${kraken_confidence} \
               --threads ${task.cpus} \
               --report ${sample_id}.merged.kraken.report \
               ${merged} \
               > ${sample_id}.merged.kraken.raw

    extract_kraken_reads.py -k ${sample_id}.merged.kraken.raw \
        --max 1000000000 --report ${sample_id}.merged.kraken.report \
        --taxid ${extract_reads_taxid} ${extract_reads_options_single} \
        --fastq-output -s ${merged} \
        -o ${sample_id}_extracted_merged.fastq

    pigz --processes ${task.cpus} ${sample_id}_extracted_merged.fastq

    # ── unmerged (now interleaved single) ───────────────────────
    ${KRAKEN2} --db ${krakendb} ${kraken_options} --confidence ${kraken_confidence} \
               --threads ${task.cpus} \
               --report ${sample_id}.unmerged.kraken.report \
               ${unmerged} \
               > ${sample_id}.unmerged.kraken.raw

    extract_kraken_reads.py -k ${sample_id}.unmerged.kraken.raw \
        --max 1000000000 --report ${sample_id}.unmerged.kraken.report \
        --taxid ${extract_reads_taxid} ${extract_reads_options_single} \
        --fastq-output -s ${unmerged} \
        -o ${sample_id}_extracted_unmerged.fastq

    pigz --processes ${task.cpus} ${sample_id}_extracted_unmerged.fastq
    """
}

process krakenresults {
    tag { }
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Results/", mode: 'copy'

    input:
        path(kraken_reports)

    output:
        path("kraken_analytic_matrix.conf_${kraken_confidence}.csv")
	path("unclassifieds_kraken_analytic_matrix.conf_${kraken_confidence}.csv")

    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o kraken_analytic_matrix.conf_${kraken_confidence}.csv
    """
}


process runbracken {
    label "microbiome"
    
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
    label "microbiome"
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
    label "microbiome"
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
    mkdir -p krona
    ktImportTaxonomy -o kraken2_taxonomy_krona.html -tax krona_db $x
    """
}
