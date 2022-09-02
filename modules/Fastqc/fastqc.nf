

process fastqc {
    tag "FASTQC on $sample_id"

    publishDir "${params.output}/fastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 

    output:
    path "${sample_id}_fastqc_logs" 


    script:
    """
    mkdir ${sample_id}_fastqc_logs
    fastqc -o ${sample_id}_fastqc_logs -f fastq -q ${reads}
    """
}


process multiqc {
    errorStrategy 'ignore'
    
    publishDir "${params.output}/multiQC", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".html") > 0) "Stats/$filename"
            else {}
        }

    
    input:
    path 'data*/*' 
    path config

    output:
    path 'multiqc_report.html'

    script:
    """
    cp $config/* .
    echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
    multiqc -v .
    """
}