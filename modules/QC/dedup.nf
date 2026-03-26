process DeduplicateReadsSeqkit {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode:'copy', pattern:'*.fastq.gz'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), path(input_fq)

    output:
        tuple val(sample_id), path("${sample_id}.dedup.fastq.gz"), emit: dedup_fq
        path("${sample_id}.dedupe_seqkit.stats.log"),              emit: dedupe_stats

    script:
    """
    seqkit rmdup --threads ${task.cpus} --by-seq \
        -o ${sample_id}.dedup.fastq.gz \
        ${input_fq} \
        > ${sample_id}.dedupe_seqkit.stats.log 2>&1
    """
}