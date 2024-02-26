

process fastqc {
    tag "FASTQC on $sample_id"
    label "fastqc"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/QC_analysis/FastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 

    output:
    path "${sample_id}_fastqc_logs" 


    script:
    """
    mkdir -p ${sample_id}_fastqc_logs
    fastqc -o ${sample_id}_fastqc_logs -f fastq -q ${reads}
    """
}


process multiqc {
    tag "Running multiQC"
    label "fastqc"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/QC_analysis/", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf("multiqc_data/*") > 0) "MultiQC_stats/multiqc_data/$filename"
            else if(filename.indexOf("general_stats.txt") > 0) "MultiQC_stats/$filename"
            else if(filename.indexOf("_report.html") > 0) "MultiQC_stats/$filename"
            else {}
        }
    
    input:
        path 'data*/*' 
        path config

    output:
        path 'AMR-Bioinformatic-pipeline_multiqc_report.html'
        path 'multiqc_general_stats.txt'
        path 'AMR-Bioinformatic-pipeline_multiqc_report_data/'

    script:
    """
    cp $config/* .
    multiqc -v data* --interactive -f --cl-config "max_table_rows: 3000"
    mv AMR-Bioinformatic-pipeline_multiqc_report_data/multiqc_general_stats.txt .
    """
}
