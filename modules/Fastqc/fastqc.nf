

process fastqc {
    tag "FASTQC on $sample_id"
    label "fastqc"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/QC_analysis/FastQC", mode: 'copy'

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
    label "fastqc"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/QC_analysis/", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf("multiqc_data/*") > 0) "MultiQC_stats/multiqc_data/$filename"
            if(filename.indexOf("general_stats.txt") > 0) "MultiQC_stats/$filename"
            if(filename.indexOf("_report.html") > 0) "MultiQC_stats/$filename"
            else {}
        }
    
    input:
    path 'data*/*' 
    path config

    output:
    path 'multiqc_report.html'
    path 'multiqc_general_stats.txt'

    script:
    """
    cp $config/* .
    multiqc -v .
    mv multiqc_data/multiqc_general_stats.txt .
    """
}
