process PE_DeduplicateReadsSeqkit {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'

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

process PE_DeduplicateReadsClumpify {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id),
              path("${sample_id}_R1.dedup.fastq.gz"),
              path("${sample_id}_R2.dedup.fastq.gz"), emit: dedup_pe_fq
        path("${sample_id}.dedupe_clumpify.stats.log"), emit: dedupe_stats

    script:
    // Reserve ~75% of allocated memory for the JVM heap, use "clumpify_mem_gb" variable 
    // or default to 8gb
    def jvm_mem = task.memory ? (task.memory.toGiga() * 0.75).intValue() 
                         : (params.clumpify_mem_gb ?: 8)
    """
    clumpify.sh \
        -Xmx${jvm_mem}g \
        threads=${task.cpus} \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=${sample_id}_R1.dedup.fastq.gz \
        out2=${sample_id}_R2.dedup.fastq.gz \
        dedupe=t \
        dupesubs=0 \
        > ${sample_id}.dedupe_clumpify.stats.log 2>&1
    """
}


process PE_DeduplicateMergedReadsClumpify {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'

    input:
        tuple val(sample_id), path(reads)   // reads[0]=merged, reads[1]=unmerged

    output:
        tuple val(sample_id),
              path("${sample_id}.merged.dedup.fastq.gz"),
              path("${sample_id}.unmerged.dedup.fastq.gz"), emit: dedup_merged_fq
        path("${sample_id}.dedupe_clumpify.stats.log"),      emit: dedupe_stats

    script:
    // Reserve ~75% of allocated memory for the JVM heap, use "clumpify_mem_gb" variable 
    // or default to 8gb
    def jvm_mem = task.memory ? (task.memory.toGiga() * 0.75).intValue() 
                         : (params.clumpify_mem_gb ?: 8)
    """
    # Merged reads: SE-like, deduplicate independently
    clumpify.sh \
        -Xmx${jvm_mem}g \
        threads=${task.cpus} \
        in=${reads[0]} \
        out=${sample_id}.merged.dedup.fastq.gz \
        dedupe=t \
        dupesubs=0 \
        >  ${sample_id}.dedupe_clumpify.stats.log 2>&1

    # Unmerged reads: also SE-like, deduplicate independently
    clumpify.sh \
        -Xmx${jvm_mem}g \
        threads=${task.cpus} \
        in=${reads[1]} \
        out=${sample_id}.unmerged.dedup.fastq.gz \
        dedupe=t \
        dupesubs=0 \
        >> ${sample_id}.dedupe_clumpify.stats.log 2>&1
    """
}


process SE_DeduplicateReadsClumpify {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'

    input:
        tuple val(sample_id), path(read)

    output:
        tuple val(sample_id),
              path("${sample_id}.dedup.fastq.gz"),          emit: dedup_fq
        path("${sample_id}.dedupe_clumpify.stats.log"),      emit: dedupe_stats

    script:
    // Reserve ~75% of allocated memory for the JVM heap, use "clumpify_mem_gb" variable 
    // or default to 8gb
    def jvm_mem = task.memory ? (task.memory.toGiga() * 0.75).intValue() 
                         : (params.clumpify_mem_gb ?: 8)
    """
    clumpify.sh \
        -Xmx${jvm_mem}g \
        threads=${task.cpus} \
        in=${read} \
        out=${sample_id}.dedup.fastq.gz \
        dedupe=t \
        dupesubs=0 \
        > ${sample_id}.dedupe_clumpify.stats.log 2>&1
    """
}