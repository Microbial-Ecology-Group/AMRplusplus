

process fastqc {
    tag "FASTQC on $sample_id"
    label "micro"

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
    label "micro"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/QC_analysis/", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf("general_stats.txt") > 0) "MultiQC_stats/$filename"
            else if(filename.indexOf("_report.html") > 0) "MultiQC_stats/$filename"
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
    multiqc -v data* --interactive -f --cl-config "max_table_rows: 5000" --outdir multiqc_data --filename multiqc_report.html
    mv multiqc_data/multiqc_report_data/multiqc_general_stats.txt .
    mv multiqc_data/multiqc_report.html .

    """
}