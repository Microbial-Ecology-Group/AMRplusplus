threads = params.threads

/*
 * Merge overlapping PE reads with FLASH
 * ------------------------------------------------------------
 *  ▸ takes one tuple  ( sample_id , [ R1.fastq.gz , R2.fastq.gz ] )
 *  ▸ runs:  flash  -M 120  -o ${sample_id}  --interleaved-output  -z  -t ${task.cpus}
 *  ▸ publishes the three FLASH FASTQ‑gz outputs + log
 */
process MergeReadsFlash {

    tag   { sample_id }
    label "small"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Flash_reads", mode: 'copy', pattern: '*.fastq.gz',
        saveAs: { fn -> fn }                       // keep original FLASH filenames

    input:
        tuple val(sample_id), path(reads)          // reads[0] = R1, reads[1] = R2

    output:
        tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), emit: merged
        tuple val(sample_id), path("${sample_id}.notCombined.fastq.gz"),  emit: unmerged
        path("${sample_id}.hist"),                                                    emit: hist
        path("${sample_id}.log"),                                                     emit: flash_log

    script:
    """
    flash -M 120 -o ${sample_id} --interleaved-output -z -t ${task.cpus} ${reads[0]} ${reads[1]}  &> ${sample_id}.log
    """
}

process SeqkitReadCounts {

    tag "seqkit"
    label "micro"

    publishDir "${params.output}/Results/Stats", mode: 'copy'

    /*
     * fastqs  – a *list* of Path objects (all merged & unmerged reads)
     * prefix  – e.g. params.QC_prefix  →  <prefix>_reads.txt
     */
    input:
        path fastqs            // list because we used collect()
        val  prefix

    output:
        path("${prefix}_reads.txt"), emit: readCounts

    script:
    """
    set -euo pipefail
    seqkit stats -T -j ${task.cpus} ${fastqs.join(' ')} \
        > ${prefix}_reads.txt
    """
}