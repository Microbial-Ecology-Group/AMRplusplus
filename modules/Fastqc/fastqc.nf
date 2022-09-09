

process fastqc {
    tag "FASTQC on $sample_id"
    conda = "$baseDir/envs/fastqc.yaml"
    container = 'enriquedoster/amrplusplus_qc:latest'

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
    conda = "$baseDir/envs/fastqc.yaml"
    container = 'enriquedoster/amrplusplus_qc:latest'
    
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
    path 'multiqc_data/multiqc_general_stats.txt'

    script:
    """
    cp $config/* .
    multiqc -v .
    """
}
