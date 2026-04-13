process SE_DeduplicateReadsSeqkit {
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

process PE_DeduplicateReadsSeqkit {
    tag   { sample_id }
    label "medium"
    publishDir "${params.output}/Deduped_reads",
               mode: 'copy', pattern: '*.fastq.gz'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id),
              path("deduped_files/${sample_id}_R1.dedup.fastq.gz"),
              path("deduped_files/${sample_id}_R2.dedup.fastq.gz"), emit: dedup_pe_fq
        path("${sample_id}.dedupe_seqkit.stats.log"),               emit: dedupe_stats

    script:
    """
    mkdir -p deduped_files

    seqkit rmdup \
        --threads ${task.cpus} \
        --by-seq \
        -1 ${r1} \
        -2 ${r2} \
        --out-dir deduped_files \
        > ${sample_id}.dedupe_seqkit.stats.log 2>&1

    # seqkit rmdup names outputs after the inputs — rename to convention
    mv deduped_files/${r1.baseName}.* deduped_files/${sample_id}_R1.dedup.fastq.gz || true
    mv deduped_files/${r2.baseName}.* deduped_files/${sample_id}_R2.dedup.fastq.gz || true
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
        tuple val(sample_id), path(merged_fq), path(unmerged_fq)

    output:
        tuple val(sample_id),
              path("${sample_id}.merged.dedup.fastq.gz"),
              path("${sample_id}.unmerged.dedup.fastq.gz"), emit: dedup_merged_fq
        path("${sample_id}.dedupe_seqkit.stats.log"),        emit: dedupe_stats

    script:
    """
    # Merged reads: fully SE in nature, deduplicate independently
    seqkit rmdup \
        --threads ${task.cpus} \
        --by-seq \
        -o ${sample_id}.merged.dedup.fastq.gz \
        ${merged_fq} \
        >  ${sample_id}.dedupe_seqkit.stats.log 2>&1

    # Unmerged reads: also SE in nature
    seqkit rmdup \
        --threads ${task.cpus} \
        --by-seq \
        -o ${sample_id}.unmerged.dedup.fastq.gz \
        ${unmerged_fq} \
        >> ${sample_id}.dedupe_seqkit.stats.log 2>&1
    """
}