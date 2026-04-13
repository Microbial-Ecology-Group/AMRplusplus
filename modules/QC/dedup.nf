process PE_DeduplicateReadsSeqkit {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id),
              path("${sample_id}_R1.dedup.fastq.gz"),
              path("${sample_id}_R2.dedup.fastq.gz"), emit: dedup_pe_fq
        path("${sample_id}.dedupe_seqkit.stats.log"),  emit: dedupe_stats

    script:
    """
    # Deduplicate R1
    seqkit rmdup \
        -j ${task.cpus} \
        -s \
        -o ${sample_id}_R1.dedup.fastq.gz \
        ${reads[0]} \
        >  ${sample_id}.dedupe_seqkit.stats.log 2>&1

    # Deduplicate R2 independently
    seqkit rmdup \
        -j ${task.cpus} \
        -s \
        -o ${sample_id}_R2.dedup.fastq.gz \
        ${reads[1]} \
        >> ${sample_id}.dedupe_seqkit.stats.log 2>&1
    """
}


process PE_DeduplicateMergedReadsSeqkit {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), path(reads)   // reads[0]=merged, reads[1]=unmerged

    output:
        tuple val(sample_id),
              path("${sample_id}.merged.dedup.fastq.gz"),
              path("${sample_id}.unmerged.dedup.fastq.gz"), emit: dedup_merged_fq
        path("${sample_id}.dedupe_seqkit.stats.log"),        emit: dedupe_stats

    script:
    """
    # Deduplicate merged reads
    seqkit rmdup \
        -j ${task.cpus} \
        -s \
        -o ${sample_id}.merged.dedup.fastq.gz \
        ${reads[0]} \
        >  ${sample_id}.dedupe_seqkit.stats.log 2>&1

    # Deduplicate unmerged reads
    seqkit rmdup \
        -j ${task.cpus} \
        -s \
        -o ${sample_id}.unmerged.dedup.fastq.gz \
        ${reads[1]} \
        >> ${sample_id}.dedupe_seqkit.stats.log 2>&1
    """
}


process SE_DeduplicateReadsSeqkit {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), path(read)

    output:
        tuple val(sample_id),
              path("${sample_id}.dedup.fastq.gz"),     emit: dedup_fq
        path("${sample_id}.dedupe_seqkit.stats.log"),  emit: dedupe_stats

    script:
    """
    seqkit rmdup \
        -j ${task.cpus} \
        -s \
        -o ${sample_id}.dedup.fastq.gz \
        ${read} \
        > ${sample_id}.dedupe_seqkit.stats.log 2>&1
    """
}