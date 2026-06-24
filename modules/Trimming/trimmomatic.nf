include {adapter_error} from './modules/nf-functions.nf'

process runqc {
    tag { sample_id }
    label "micro_long"

    publishDir "${params.output}/QC_trimming", mode: 'copy', pattern: '*.fastq.gz',
        saveAs: { filename ->
            if(filename.indexOf("P.fastq.gz") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq.gz") > 0) "Unpaired/$filename"
            else {}
        }

    input:
        tuple val(sample_id), path(reads)  

    output:
        tuple val(sample_id), path("${sample_id}*P.fastq.gz"), emit: paired_fastq
        tuple val(sample_id), path("${sample_id}*U.fastq.gz"), emit: unpaired_fastq
        path("${sample_id}.trimmomatic.stats.log"), emit: trimmomatic_stats
   script:
    """
     ${TRIMMOMATIC} \
      PE \
      -threads ${task.cpus} \
      ${reads[0]} ${reads[1]} ${sample_id}.1P.fastq.gz ${sample_id}.1U.fastq.gz ${sample_id}.2P.fastq.gz ${sample_id}.2U.fastq.gz \
      ILLUMINACLIP:${params.adapters}:2:30:10:3:TRUE \
      LEADING:${params.leading} \
      TRAILING:${params.trailing} \
      SLIDINGWINDOW:${params.slidingwindow} \
      MINLEN:${params.minlen} \
      CROP:${params.crop_len} \
      2> ${sample_id}.trimmomatic.stats.log
      
    """
}

process runqc_se {
  tag { sample_id }
  label "small"

  publishDir "${params.output}/QC_trimming/Single", mode: 'copy', pattern: '*.fastq.gz',
    saveAs: { fn -> fn }

  input:
    tuple val(sample_id), path(read)

  output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"),            emit: se_fastq
    path("${sample_id}.trimmomatic.stats.log"),                              emit: trimmomatic_stats    // keep if you still want the stderr log
    path("${sample_id}.trimmomatic.summary.txt"),                            emit: trimmomatic_summary  // NEW: uniform summary

  script:
  """
  ${TRIMMOMATIC} \
    SE \
    -threads ${task.cpus} \
    -summary ${sample_id}.trimmomatic.summary.txt \
    ${read} ${sample_id}.trimmed.fastq.gz \
    ILLUMINACLIP:${params.adapters}:2:30:10:3:TRUE \
    LEADING:${params.leading} \
    TRAILING:${params.trailing} \
    SLIDINGWINDOW:${params.slidingwindow} \
    MINLEN:${params.minlen} \
    CROP:${params.crop_len} \
    2> ${sample_id}.trimmomatic.stats.log
  """
}

process QCstats {
    tag "Make QC summary file"
    label "small"

    publishDir "${params.output}/Results", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
            else {}
        }

    input:
        file(stats)

    output:
        path("trimmomatic.stats"), emit: combo_trim_stats

   script:
    """
    ${PYTHON3} $baseDir/bin/trimmomatic_stats.py -i ${stats} -o trimmomatic.stats
    """
}

process QCstats_SE {
    tag "Make QC summary file (SE)"
    label "small"

    publishDir "${params.output}/Results", mode: 'copy',
       saveAs: { fn -> fn.endsWith(".stats") ? "Stats/$fn" : null }

    input:
      file(summaries)

    output:
      path("trimmomatic.stats"), emit: combo_trim_stats

    script:
    """

    """
}

