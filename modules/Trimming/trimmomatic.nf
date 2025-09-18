include {adapter_error} from "$baseDir/modules/nf-functions.nf"

if( params.adapters ) {
    adapters = file(params.adapters)
    if( !adapters.exists() ) return adapter_error(adapters)
}

threads = params.threads
min = params.min
max = params.max
skip = params.skip
samples = params.samples

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

process runqc {
    tag { sample_id }
    label "trimming"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

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

    """
     ${TRIMMOMATIC} \
      PE \
      -threads ${threads} \
      ${reads[0]} ${reads[1]} ${sample_id}.1P.fastq.gz ${sample_id}.1U.fastq.gz ${sample_id}.2P.fastq.gz ${sample_id}.2U.fastq.gz \
      ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:${slidingwindow} \
      MINLEN:${minlen} \
      2> ${sample_id}.trimmomatic.stats.log
      
    """
}

process runqc_se {
  tag { sample_id }
  label "trimming"

  publishDir "${params.output}/QC_trimming_SE", mode: 'copy', pattern: '*.fastq.gz',
    saveAs: { fn -> fn }

  input:
    tuple val(sample_id), path(read)

  output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"),            emit: se_fastq
    path("${sample_id}.trimmomatic.stats.log"),                              emit: trimmomatic_stats    // keep if you still want the stderr log
    path("${sample_id}.trimmomatic.summary.txt"),                            emit: trimmomatic_summary  // NEW: uniform summary

  """
  ${TRIMMOMATIC} \
    SE \
    -threads ${threads} \
    -summary ${sample_id}.trimmomatic.summary.txt \
    ${read} ${sample_id}.trimmed.fastq.gz \
    ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
    LEADING:${leading} \
    TRAILING:${trailing} \
    SLIDINGWINDOW:${slidingwindow} \
    MINLEN:${minlen} \
    2> ${sample_id}.trimmomatic.stats.log
  """
}

process QCstats {
    tag "Make QC summary file"
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Results", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
            else {}
        }

    input:
        file(stats)

    output:
        path("trimmomatic.stats"), emit: combo_trim_stats

    """
    ${PYTHON3} $baseDir/bin/trimmomatic_stats.py -i ${stats} -o trimmomatic.stats
    """
}

process QCstats_SE {
  tag "Make QC summary file (SE)"
  label "python"

  publishDir "${params.output}/Results", mode: 'copy',
    saveAs: { fn -> fn.endsWith(".stats") ? "Stats/$fn" : null }

  input:
    file(summaries)  // collected *.trimmomatic.summary.txt files

  output:
    path("trimmomatic.stats"), emit: combo_trim_stats

  """
    set -euo pipefail
    echo -e "sample\\ttotal\\tforward_surviving\\treverse_surviving\\tdropped" > trimmomatic.stats
    for f in ${summaries}; do
      sample=`basename "\\$f" .trimmomatic.summary.txt`
      line=`grep -E '^Input Reads' "\\$f" | tail -1 || true`
      total=`echo "\\$line" | grep -oE 'Input Reads: *[0-9,]+' | awk '{gsub(/,/,"",\\$3); print \\$3}'`
      surv=`echo "\\$line"  | grep -oE 'Surviving: *[0-9,]+'   | awk '{gsub(/,/,"",\\$2); print \\$2}'`
      drop=`echo "\\$line"  | grep -oE 'Dropped: *[0-9,]+'     | awk '{gsub(/,/,"",\\$2); print \\$2}'`
      total=\\${total:-0}; surv=\\${surv:-0}; drop=\\${drop:-0}
      echo -e "\\$sample\\t\\$total\\t\\$surv\\t0\\t\\$drop" >> trimmomatic.stats
    done
  """
}
