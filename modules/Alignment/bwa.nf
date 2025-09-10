include { reference_error ; amr_error ; annotation_error } from "$baseDir/modules/nf-functions.nf"


if( params.amr ) {
    amr = file(params.amr)
    if( !amr.exists() ) return amr_error(amr)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}


threads = params.threads

deduped = params.deduped

process index {
    tag "Creating bwa index"
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/BWA_Index", mode: "copy"

    input:
    path fasta

    output: 
    path("${fasta}*"), emit: bwaindex, includeInputs: true

    script:
    """
    bwa index ${fasta}
    #--threads $task.cpus 
    """
}


process bwa_align {
    tag "$pair_id"
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/BAM_files", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("_alignment_sorted.bam") > 0) "Standard/$filename"
            else if(filename.indexOf("_alignment_dedup.bam") > 0) "Deduped/$filename"
            else {}
        }

    input:
        path indexfiles 
        tuple val(pair_id), path(reads) 

    output:
        tuple val(pair_id), path("${pair_id}_alignment_sorted.bam"), emit: bwa_bam
        tuple val(pair_id), path("${pair_id}_alignment_dedup.bam"), emit: bwa_dedup_bam, optional: true

    script:
    if( deduped == "N")
        """
        ${BWA} mem ${indexfiles[0]} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
        ${SAMTOOLS} view -@ ${threads} -S -b ${pair_id}_alignment.sam > ${pair_id}_alignment.bam
        rm ${pair_id}_alignment.sam
        ${SAMTOOLS} sort -@ ${threads} -n ${pair_id}_alignment.bam -o ${pair_id}_alignment_sorted.bam
        rm ${pair_id}_alignment.bam
        """
    else if( deduped == "Y")
        """
        ${BWA} mem ${indexfiles[0]} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
        ${SAMTOOLS} view -@ ${threads} -S -b ${pair_id}_alignment.sam > ${pair_id}_alignment.bam
        rm ${pair_id}_alignment.sam
        ${SAMTOOLS} sort -@ ${threads} -n ${pair_id}_alignment.bam -o ${pair_id}_alignment_sorted.bam
        rm ${pair_id}_alignment.bam
        ${SAMTOOLS} fixmate -@ ${threads} ${pair_id}_alignment_sorted.bam ${pair_id}_alignment_sorted_fix.bam
        ${SAMTOOLS} sort -@ ${threads} ${pair_id}_alignment_sorted_fix.bam -o ${pair_id}_alignment_sorted_fix.sorted.bam
        rm ${pair_id}_alignment_sorted_fix.bam
        ${SAMTOOLS} rmdup -S ${pair_id}_alignment_sorted_fix.sorted.bam ${pair_id}_alignment_dedup.bam
        rm ${pair_id}_alignment_sorted_fix.sorted.bam
        ${SAMTOOLS} view -@ ${threads} -h -o ${pair_id}_alignment_dedup.sam ${pair_id}_alignment_dedup.bam
        rm ${pair_id}_alignment_dedup.sam
        """
    else
        error "Invalid deduplication flag --deduped: ${deduped}. Please use --deduped Y for deduplicated counts, or avoid using this flag altogether to skip this error."
}

process bwa_merged_align {

    tag   { sample_id }
    label "alignment"

    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    publishDir "${params.output}/Alignment/BAM_files", mode: 'copy',
        saveAs: { fn ->
            // “Standard” = sorted-only BAMs, “Deduped” = after rmdup
            if (fn.endsWith('_merged_alignment_sorted.bam')   ||
                fn.endsWith('_unmerged_alignment_sorted.bam'))
                    "Standard/$fn"
            else if (fn.endsWith('_merged_alignment_dedup.bam') ||
                     fn.endsWith('_unmerged_alignment_dedup.bam'))
                    "Deduped/$fn"
            else
                    null
        }

    /* ───────── inputs ───────────────────────────────────────────────
     *  indexfiles[0]  – BWA index prefix (6 files)
     *  merged_fq      – FLASH-merged single-end reads
     *  unmerged_fq    – FLASH-unmerged single-end reads
     */
    input:
        path  indexfiles
        tuple val(sample_id), path(merged_fq), path(unmerged_fq)

    /* ───────── outputs ────────────────────────────────────────────── */
    output:
        tuple val(sample_id), path("${sample_id}_merged_alignment_sorted.bam"),   emit: merged_bam
        tuple val(sample_id), path("${sample_id}_unmerged_alignment_sorted.bam"), emit: unmerged_bam

        tuple val(sample_id), path("${sample_id}_merged_alignment_dedup.bam"),    emit: merged_dedup_bam,    optional: true
        tuple val(sample_id), path("${sample_id}_unmerged_alignment_dedup.bam"),  emit: unmerged_dedup_bam,  optional: true

    /* ───────── script ─────────────────────────────────────────────── */
    script:
    def cpu = task.cpus ?: threads                            // convenience
    if (deduped == 'N') """
        set -euo pipefail

        # ───── merged reads ──────────────────────────────────────────
        bwa mem ${indexfiles[0]} ${merged_fq} -t ${cpu} \
             -R '@RG\\tID:${sample_id}_merged\\tSM:${sample_id}' \
        | samtools sort -@ ${cpu} -n -o ${sample_id}_merged_alignment_sorted.bam -

        # ───── un-merged reads ───────────────────────────────────────
        bwa mem ${indexfiles[0]} ${unmerged_fq} -t ${cpu} \
             -R '@RG\\tID:${sample_id}_unmerged\\tSM:${sample_id}' \
        | samtools sort -@ ${cpu} -n -o ${sample_id}_unmerged_alignment_sorted.bam -
    """
    else if (deduped == 'Y') """
        set -euo pipefail

        # ───── merged reads (deduped) ────────────────────────────────
        bwa mem ${indexfiles[0]} ${merged_fq} -t ${cpu} \
             -R '@RG\\tID:${sample_id}_merged\\tSM:${sample_id}' \
        | samtools sort -@ ${cpu} -n -o ${sample_id}_merged_alignment_sorted.bam -

        samtools fixmate -@ ${cpu} ${sample_id}_merged_alignment_sorted.bam tmp.bam
        samtools sort    -@ ${cpu} tmp.bam -o tmp.srt.bam
        samtools rmdup -S tmp.srt.bam ${sample_id}_merged_alignment_dedup.bam
        rm tmp.bam tmp.srt.bam

        # ───── un-merged reads (deduped) ─────────────────────────────
        bwa mem ${indexfiles[0]} ${unmerged_fq} -t ${cpu} \
             -R '@RG\\tID:${sample_id}_unmerged\\tSM:${sample_id}' \
        | samtools sort -@ ${cpu} -n -o ${sample_id}_unmerged_alignment_sorted.bam -

        samtools fixmate -@ ${cpu} ${sample_id}_unmerged_alignment_sorted.bam tmp.bam
        samtools sort    -@ ${cpu} tmp.bam -o tmp.srt.bam
        samtools rmdup -S tmp.srt.bam ${sample_id}_unmerged_alignment_dedup.bam
        rm tmp.bam tmp.srt.bam
    """
    else
        error "Invalid --deduped flag: ${deduped}. Use Y or N."
}

process bwa_align_se {
    tag { sample_id }
    label "alignment"

    publishDir "${params.output}/Alignment/BAM_files", mode: "copy",
        saveAs: { fn -> fn.endsWith("_alignment_sorted.bam") ? "Standard/$fn" : null }

    input:
        path indexfiles
        tuple val(sample_id), path(read)

    output:
        tuple val(sample_id), path("${sample_id}_alignment_sorted.bam"), emit: bwa_bam

    script:
    """
    set -euo pipefail

    ${BWA} mem ${indexfiles[0]} ${read} -t ${task.cpus} \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' \
      | ${SAMTOOLS} sort -@ ${task.cpus} -o ${sample_id}_alignment_sorted.bam -

    ${SAMTOOLS} index ${sample_id}_alignment_sorted.bam
    """
}

process samtools_dedup_se {
    tag { sample_id }
    label "alignment"

    publishDir "${params.output}/Alignment/BAM_files/Deduped", mode: "copy"

    input:
        tuple val(sample_id), path(bam_in)  // coord-sorted BAM with index

    output:
        tuple val(sample_id), path("${sample_id}_alignment_dedup.bam"), emit: dedup_bam

    script:
    """
    set -euo pipefail

    # ensure coordinate sort (safe to re-sort; markdup requires coord-sorted)
    ${SAMTOOLS} sort -@ ${task.cpus} -o ${sample_id}.coord.bam ${bam_in}
    ${SAMTOOLS} index ${sample_id}.coord.bam

    # remove duplicates
    ${SAMTOOLS} markdup -r -@ ${task.cpus} ${sample_id}.coord.bam ${sample_id}_alignment_dedup.bam
    ${SAMTOOLS} index ${sample_id}_alignment_dedup.bam

    # clean temp
    rm -f ${sample_id}.coord.bam ${sample_id}.coord.bam.bai
    """
}


process bwa_rm_contaminant_fq {
    tag { pair_id }
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3 
 
    publishDir "${params.output}/HostRemoval", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("fastq.gz") > 0) "NonHostFastq/$filename"
            else {}
        }

    input:
    path indexfiles
    tuple val(pair_id), path(reads) 

    output:
    tuple val(pair_id), path("${pair_id}.non.host.R*.fastq.gz"), emit: nonhost_reads
    path("${pair_id}.samtools.idxstats"), emit: host_rm_stats
    
    """
    ${BWA} mem ${indexfiles[0]} ${reads[0]} ${reads[1]} -t ${threads} > ${pair_id}.host.sam
    ${SAMTOOLS} view -bS ${pair_id}.host.sam | ${SAMTOOLS} sort -@ ${threads} -o ${pair_id}.host.sorted.bam
    rm ${pair_id}.host.sam
    ${SAMTOOLS} index ${pair_id}.host.sorted.bam && ${SAMTOOLS} idxstats ${pair_id}.host.sorted.bam > ${pair_id}.samtools.idxstats
    ${SAMTOOLS} view -h -f 12 -b ${pair_id}.host.sorted.bam -o ${pair_id}.host.sorted.removed.bam
    ${SAMTOOLS} sort -n -@ ${threads} ${pair_id}.host.sorted.removed.bam -o ${pair_id}.host.resorted.removed.bam
    ${SAMTOOLS}  \
       fastq -@ ${threads} -c 6  \
      ${pair_id}.host.resorted.removed.bam \
      -1 ${pair_id}.non.host.R1.fastq.gz \
      -2 ${pair_id}.non.host.R2.fastq.gz \
      -0 /dev/null -s /dev/null -n

    rm *.bam
    """

}

process bwa_rm_contaminant_merged_fq {

    tag   { sample_id }
    label "medium"

    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    publishDir "${params.output}/HostRemoval", mode: 'copy',
        saveAs: { fn ->
            if (fn.endsWith('.fastq.gz'))
                "NonHostFastq/$fn"          // keep only the de-contaminated FASTQs
            else
                null
        }

    /* ───────── inputs ────────────────────────────────────────────────
     *  indexfiles[0]  – bwa index prefix (6 files in the same dir)
     *  merged_fq      – FLASH-merged reads      (sample.extendedFrags.fastq.gz)
     *  unmerged_fq    – FLASH-unmerged reads    (sample.notCombined.fastq.gz)
     */
    input:
        path  indexfiles
        tuple val(sample_id), path(merged_fq), path(unmerged_fq)

    /* ───────── outputs ─────────────────────────────────────────────── */
    output:
        tuple val(sample_id), path("${sample_id}.merged.non.host.fastq.gz"),   emit: nonhost_merged
        path("${sample_id}.merged.samtools.idxstats"),                         emit: host_rm_stats_merged

        tuple val(sample_id), path("${sample_id}.unmerged.non.host.fastq.gz"), emit: nonhost_unmerged
        path("${sample_id}.unmerged.samtools.idxstats"),                       emit: host_rm_stats_unmerged

    script:
    """
    set -euo pipefail

    # ───────────────────────── merged reads ────────────────────────────
    bwa mem ${indexfiles[0]} ${merged_fq} -t ${task.cpus} \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.merged.host.sorted.bam

    samtools index   ${sample_id}.merged.host.sorted.bam
    samtools idxstats ${sample_id}.merged.host.sorted.bam \
        > ${sample_id}.merged.samtools.idxstats

    samtools view -b -f 4 ${sample_id}.merged.host.sorted.bam \
      | samtools fastq -@ ${task.cpus} -c 6 - | pigz -p ${task.cpus} -c > ${sample_id}.merged.non.host.fastq.gz


    # ──────────────────────── un-merged reads ──────────────────────────
    bwa mem ${indexfiles[0]} ${unmerged_fq} -t ${task.cpus} \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.unmerged.host.sorted.bam

    samtools index   ${sample_id}.unmerged.host.sorted.bam
    samtools idxstats ${sample_id}.unmerged.host.sorted.bam \
        > ${sample_id}.unmerged.samtools.idxstats

    samtools view -b -f 4 ${sample_id}.unmerged.host.sorted.bam \
      | samtools fastq -@ ${task.cpus} -c 6 - | pigz -p ${task.cpus} -c > ${sample_id}.unmerged.non.host.fastq.gz
    """
}

process bwa_rm_contaminant_se {
    tag { sample_id }
    label "alignment"

    publishDir "${params.output}/HostRemoval", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.fastq.gz') ? "NonHostFastq/$fn" : null }

    input:
        path  indexfiles
        tuple val(sample_id), path(read)

    output:
        tuple val(sample_id), path("${sample_id}.non.host.fastq.gz"), emit: nonhost_reads
        path("${sample_id}.samtools.idxstats"),                        emit: host_rm_stats

    script:
    """
    set -euo pipefail

    bwa mem ${indexfiles[0]} ${read} -t ${task.cpus} \
      | samtools sort -@ ${task.cpus} -o ${sample_id}.host.sorted.bam -

    samtools index ${sample_id}.host.sorted.bam
    samtools idxstats ${sample_id}.host.sorted.bam > ${sample_id}.samtools.idxstats

    # keep only reads UNMAPPED to host
    samtools view -b -f 4 ${sample_id}.host.sorted.bam \
      | samtools fastq -@ ${task.cpus} -c 6 - \
      | pigz -p ${task.cpus} -c > ${sample_id}.non.host.fastq.gz
    """
}

process HostRemovalStats {
    tag { sample_id }
    label "alignment"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3 

    publishDir "${params.output}/Results", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
        }

    input:
        file(host_rm_stats)

    output:
        path("host.removal.stats"), emit: combo_host_rm_stats

    """
    ${PYTHON3} $baseDir/bin/samtools_idxstats.py -i ${host_rm_stats} -o host.removal.stats
    """
}

process samtools_merge_bams {
    tag { sample_id }
    label 'alignment'

    publishDir "${params.output}/Alignment/BAM_files/Combined", mode: 'copy'

    input:
        tuple val(sample_id), path(bam_list)  // list‑of‑two BAMs

    output:
        tuple val(sample_id), path("${sample_id}_combined.bam"), emit: combo_bam

    script:
    def cpu = task.cpus ?: 4
    """
    samtools merge -@ ${cpu} ${sample_id}_combined.unsorted.bam ${bam_list.join(' ')}

    samtools sort  -@ ${cpu} -o ${sample_id}_combined.bam  ${sample_id}_combined.unsorted.bam

    samtools index ${sample_id}_combined.bam
    """
}