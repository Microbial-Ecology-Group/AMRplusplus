process index {
    tag "Creating bwa index"
    label "micro"

    publishDir "${params.output}/Alignment/BWA_Index", mode: "copy"

    input:
    path fasta

    output: 
    path("${fasta}*"), emit: bwaindex, includeInputs: true

    script:
    """
    \$BWA index ${fasta}
    #--threads $task.cpus 
    """
}


process bwa_align {
    tag "$pair_id"
    label "small"

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
    if( params.deduped == "N")
        """
        \$BWA mem ${indexfiles[0]} ${reads} -t ${task.cpus} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
        \$SAMTOOLS view -@ ${task.cpus} -S -b ${params.samtools_flag} ${pair_id}_alignment.sam > ${pair_id}_alignment.bam
        rm ${pair_id}_alignment.sam
        \$SAMTOOLS sort -@ ${task.cpus} -n ${pair_id}_alignment.bam -o ${pair_id}_alignment_sorted.bam
        rm ${pair_id}_alignment.bam
        """
    else if( params.deduped == "Y")
        """
        \$BWA mem ${indexfiles[0]} ${reads} -t ${task.cpus} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}_alignment.sam
        \$SAMTOOLS view -@ ${task.cpus} -S -b ${params.samtools_flag} ${pair_id}_alignment.sam > ${pair_id}_alignment.bam
        rm ${pair_id}_alignment.sam
        \$SAMTOOLS sort -@ ${task.cpus} -n ${pair_id}_alignment.bam -o ${pair_id}_alignment_sorted.bam
        rm ${pair_id}_alignment.bam
        \$SAMTOOLS fixmate -@ ${task.cpus} -m ${pair_id}_alignment_sorted.bam ${pair_id}_alignment_sorted_fix.bam
        \$SAMTOOLS sort -@ ${task.cpus} ${pair_id}_alignment_sorted_fix.bam -o ${pair_id}_alignment_sorted_fix.sorted.bam
        rm ${pair_id}_alignment_sorted_fix.bam
        \$SAMTOOLS markdup -r ${pair_id}_alignment_sorted_fix.sorted.bam ${pair_id}_alignment_dedup.bam
        rm ${pair_id}_alignment_sorted_fix.sorted.bam
        """
    else
        error "Invalid deduplication flag --deduped: ${params.deduped}. Please use --deduped Y for deduplicated counts, or avoid using this flag altogether to skip this error."
}

process bwa_merged_align {
    tag   { sample_id }
    label "small"

    publishDir "${params.output}/Alignment/BAM_files", mode: 'copy',
        saveAs: { fn ->
            if (fn.endsWith('_merged_alignment_sorted.bam')   ||
                fn.endsWith('_unmerged_alignment_sorted.bam'))
                    "Standard/$fn"
            else if (fn.endsWith('_merged_alignment_dedup.bam') ||
                     fn.endsWith('_unmerged_alignment_dedup.bam'))
                    "Deduped/$fn"
            else
                    null
        }

    input:
        path  indexfiles
        tuple val(sample_id), path(merged_fq), path(unmerged_fq)

    output:
        tuple val(sample_id), path("${sample_id}_merged_alignment_sorted.bam"),   emit: merged_bam
        tuple val(sample_id), path("${sample_id}_unmerged_alignment_sorted.bam"), emit: unmerged_bam

        tuple val(sample_id), path("${sample_id}_merged_alignment_dedup.bam"),    emit: merged_dedup_bam,    optional: true
        tuple val(sample_id), path("${sample_id}_unmerged_alignment_dedup.bam"),  emit: unmerged_dedup_bam,  optional: true

    script:
  
    if (params.deduped == 'N') """
        set -euo pipefail

        # ── helpers ──────────────────────────────────────────────────
        has_reads() { [ "\$(zcat "\$1" 2>/dev/null | head -c 1 | wc -c)" -gt 0 ]; }
        empty_bam() { printf '@HD\\tVN:1.6\\tSO:unsorted\\n' | \$SAMTOOLS view -bS -o "\$1" -; }

        # ───── merged reads ──────────────────────────────────────────
        if has_reads ${merged_fq}; then
            \$BWA mem ${indexfiles[0]} ${merged_fq} -t ${task.cpus} \
                -R '@RG\\tID:${sample_id}_merged\\tSM:${sample_id}' \
            | \$SAMTOOLS view -@ ${task.cpus} -b ${params.samtools_flag} - \
            | \$SAMTOOLS sort -@ ${task.cpus} -n -o ${sample_id}_merged_alignment_sorted.bam -
        else
            echo "[INFO] No merged reads for ${sample_id} — creating empty BAM"
            empty_bam ${sample_id}_merged_alignment_sorted.bam
        fi

        # ───── un-merged reads ───────────────────────────────────────
        if has_reads ${unmerged_fq}; then
            \$BWA mem ${indexfiles[0]} ${unmerged_fq} -t ${task.cpus} \
                -R '@RG\\tID:${sample_id}_unmerged\\tSM:${sample_id}' \
            | \$SAMTOOLS view -@ ${task.cpus} -b ${params.samtools_flag} - \
            | \$SAMTOOLS sort -@ ${task.cpus} -n -o ${sample_id}_unmerged_alignment_sorted.bam -
        else
            echo "[INFO] No unmerged reads for ${sample_id} — creating empty BAM"
            empty_bam ${sample_id}_unmerged_alignment_sorted.bam
        fi
    """
    else if (params.deduped == 'Y') """
        set -euo pipefail

        # ── helpers ──────────────────────────────────────────────────
        has_reads() { [ "\$(zcat "\$1" 2>/dev/null | head -c 1 | wc -c)" -gt 0 ]; }
        empty_bam() { printf '@HD\\tVN:1.6\\tSO:unsorted\\n' | \$SAMTOOLS view -bS -o "\$1" -; }

        # ───── merged reads (+ dedup) ────────────────────────────────
        if has_reads ${merged_fq}; then
            \$BWA mem ${indexfiles[0]} ${merged_fq} -t ${task.cpus} \
                -R '@RG\\tID:${sample_id}_merged\\tSM:${sample_id}' \
            | \$SAMTOOLS view -@ ${task.cpus} -b ${params.samtools_flag} - \
            | \$SAMTOOLS sort -@ ${task.cpus} -n -o ${sample_id}_merged_alignment_sorted.bam -

            \$SAMTOOLS fixmate -@ ${task.cpus} -m ${sample_id}_merged_alignment_sorted.bam tmp_merged.bam
            \$SAMTOOLS sort    -@ ${task.cpus} tmp_merged.bam -o tmp_merged.srt.bam
            \$SAMTOOLS markdup -r -@ ${task.cpus} tmp_merged.srt.bam ${sample_id}_merged_alignment_dedup.bam
            rm -f tmp_merged.bam tmp_merged.srt.bam
        else
            echo "[INFO] No merged reads for ${sample_id} — creating empty BAMs"
            empty_bam ${sample_id}_merged_alignment_sorted.bam
            empty_bam ${sample_id}_merged_alignment_dedup.bam
        fi

        # ───── un-merged reads (+ dedup) ─────────────────────────────
        if has_reads ${unmerged_fq}; then
            \$BWA mem ${indexfiles[0]} ${unmerged_fq} -t ${task.cpus} \
                -R '@RG\\tID:${sample_id}_unmerged\\tSM:${sample_id}' \
            | \$SAMTOOLS view -@ ${task.cpus} -b ${params.samtools_flag} - \
            | \$SAMTOOLS sort -@ ${task.cpus} -n -o ${sample_id}_unmerged_alignment_sorted.bam -

            \$SAMTOOLS fixmate -@ ${task.cpus} -m ${sample_id}_unmerged_alignment_sorted.bam tmp_unmerged.bam
            \$SAMTOOLS sort    -@ ${task.cpus} tmp_unmerged.bam -o tmp_unmerged.srt.bam
            \$SAMTOOLS markdup -r -@ ${task.cpus} tmp_unmerged.srt.bam ${sample_id}_unmerged_alignment_dedup.bam
            rm -f tmp_unmerged.bam tmp_unmerged.srt.bam
        else
            echo "[INFO] No unmerged reads for ${sample_id} — creating empty BAMs"
            empty_bam ${sample_id}_unmerged_alignment_sorted.bam
            empty_bam ${sample_id}_unmerged_alignment_dedup.bam
        fi
    """
    else
        error "Invalid --deduped flag: ${params.deduped}. Use Y or N."
}




process bwa_align_se {
    tag { sample_id }
    label "small"

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

    \$BWA mem ${indexfiles[0]} ${read} -t ${task.cpus} \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' \
    | \$SAMTOOLS view -@ ${task.cpus} -b ${params.samtools_flag} - \
    | \$SAMTOOLS sort -@ ${task.cpus} -o ${sample_id}_alignment_sorted.bam -

    \$SAMTOOLS index ${sample_id}_alignment_sorted.bam
    """
}


process samtools_dedup_se {
    tag { sample_id }
    label "small"

    publishDir "${params.output}/Alignment/BAM_files/Deduped", mode: "copy"

    input:
        tuple val(sample_id), path(bam_in)  // coord-sorted BAM with index

    output:
        tuple val(sample_id), path("${sample_id}_alignment_dedup.bam"), emit: dedup_bam

    script:
    """
    set -euo pipefail

    # ensure coordinate sort (safe to re-sort; markdup requires coord-sorted)
    \$SAMTOOLS sort -@ ${task.cpus} -o ${sample_id}.coord.bam ${bam_in}
    \$SAMTOOLS index ${sample_id}.coord.bam

    # remove duplicates
    \$SAMTOOLS markdup -r -@ ${task.cpus} ${sample_id}.coord.bam ${sample_id}_alignment_dedup.bam
    \$SAMTOOLS index ${sample_id}_alignment_dedup.bam

    # clean temp
    rm -f ${sample_id}.coord.bam ${sample_id}.coord.bam.bai
    """
}


process bwa_rm_contaminant_fq {
    tag { pair_id }
    label "medium"
 
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
    
    script:
    """
    \$BWA mem ${indexfiles[0]} ${reads[0]} ${reads[1]} -t ${task.cpus} > ${pair_id}.host.sam
    \$SAMTOOLS view -bS ${pair_id}.host.sam | \$SAMTOOLS sort -@ ${task.cpus} -o ${pair_id}.host.sorted.bam
    rm ${pair_id}.host.sam
    \$SAMTOOLS index ${pair_id}.host.sorted.bam && \$SAMTOOLS idxstats ${pair_id}.host.sorted.bam > ${pair_id}.samtools.idxstats
    \$SAMTOOLS view -h -f 12 -b ${pair_id}.host.sorted.bam -o ${pair_id}.host.sorted.removed.bam
    \$SAMTOOLS sort -n -@ ${task.cpus} ${pair_id}.host.sorted.removed.bam -o ${pair_id}.host.resorted.removed.bam
    \$SAMTOOLS  \
       fastq -@ ${task.cpus} -c 6  \
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

    publishDir "${params.output}/HostRemoval", mode: 'copy',
        saveAs: { fn ->
            if (fn.endsWith('.fastq.gz'))
                "NonHostFastq/$fn"
            else
                null
        }

    input:
        path  indexfiles
        tuple val(sample_id), path(merged_fq), path(unmerged_fq)

    output:
        tuple val(sample_id), path("${sample_id}.merged.non.host.fastq.gz"),   emit: nonhost_merged
        path("${sample_id}.merged.samtools.idxstats"),                         emit: host_rm_stats_merged

        tuple val(sample_id), path("${sample_id}.unmerged.non.host.fastq.gz"), emit: nonhost_unmerged
        path("${sample_id}.unmerged.samtools.idxstats"),                       emit: host_rm_stats_unmerged

    script:
    """
    set -euo pipefail

    # ── helper: true if the gzipped FASTQ contains at least one read ──
    has_reads() { [ "\$(zcat "\$1" 2>/dev/null | head -c 1 | wc -c)" -gt 0 ]; }

    # ───────────────────────── merged reads ────────────────────────────
    if has_reads ${merged_fq}; then
        \$BWA mem ${indexfiles[0]} ${merged_fq} -t ${task.cpus} \\
            | \$SAMTOOLS sort -@ ${task.cpus} -o ${sample_id}.merged.host.sorted.bam

        \$SAMTOOLS index   ${sample_id}.merged.host.sorted.bam
        \$SAMTOOLS idxstats ${sample_id}.merged.host.sorted.bam \\
            > ${sample_id}.merged.samtools.idxstats

        \$SAMTOOLS view -b -f 4 ${sample_id}.merged.host.sorted.bam \\
          | \$SAMTOOLS fastq -@ ${task.cpus} -c 6 - \\
          | pigz -p ${task.cpus} -c > ${sample_id}.merged.non.host.fastq.gz
    else
        echo "[INFO] No merged reads for ${sample_id} — writing empty outputs"
        echo -n | gzip > ${sample_id}.merged.non.host.fastq.gz
        printf "*\\t0\\t0\\t0\\n" > ${sample_id}.merged.samtools.idxstats
    fi

    # ──────────────────────── un-merged reads ──────────────────────────
    if has_reads ${unmerged_fq}; then
        \$BWA mem ${indexfiles[0]} ${unmerged_fq} -t ${task.cpus} \\
            | \$SAMTOOLS sort -@ ${task.cpus} -o ${sample_id}.unmerged.host.sorted.bam

        \$SAMTOOLS index   ${sample_id}.unmerged.host.sorted.bam
        \$SAMTOOLS idxstats ${sample_id}.unmerged.host.sorted.bam \\
            > ${sample_id}.unmerged.samtools.idxstats

        \$SAMTOOLS view -b -f 4 ${sample_id}.unmerged.host.sorted.bam \\
          | \$SAMTOOLS fastq -@ ${task.cpus} -c 6 - \\
          | pigz -p ${task.cpus} -c > ${sample_id}.unmerged.non.host.fastq.gz
    else
        echo "[INFO] No unmerged reads for ${sample_id} — writing empty outputs"
        echo -n | gzip > ${sample_id}.unmerged.non.host.fastq.gz
        printf "*\\t0\\t0\\t0\\n" > ${sample_id}.unmerged.samtools.idxstats
    fi
    """
}


process bwa_rm_contaminant_se {
    tag { sample_id }
    label "medium"

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

    \$BWA mem ${indexfiles[0]} ${read} -t ${task.cpus} \
      | \$SAMTOOLS sort -@ ${task.cpus} -o ${sample_id}.host.sorted.bam -

    \$SAMTOOLS index ${sample_id}.host.sorted.bam
    \$SAMTOOLS idxstats ${sample_id}.host.sorted.bam > ${sample_id}.samtools.idxstats

    # keep only reads UNMAPPED to host
    \$SAMTOOLS view -b -f 4 ${sample_id}.host.sorted.bam \
      | \$SAMTOOLS fastq -@ ${task.cpus} -c 6 - \
      | pigz -p ${task.cpus} -c > ${sample_id}.non.host.fastq.gz
    """
}

process HostRemovalStats {
    tag "Summarizing host removal stats"
    label "micro"

    publishDir "${params.output}/Results", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
        }

    input:
        file(host_rm_stats)

    output:
        path("host.removal.stats"), emit: combo_host_rm_stats

    script:
    """
    \$PYTHON3 $baseDir/bin/samtools_idxstats.py -i ${host_rm_stats} -o host.removal.stats
    """
}

process samtools_merge_bams {
    tag { sample_id }
    label 'small'

    publishDir { "${params.output}/Alignment/BAM_files/Combined/${subdir}" }, mode: 'copy'

    input:
        tuple val(sample_id), path(bam_list)
        val(subdir)                          // "" standard, "Deduped" for dedup
    output:
        tuple val(sample_id), path("${sample_id}_combined*.bam*"), emit: combo_bam
    script:
    def outname = subdir ? "${sample_id}_combined.dedup" : "${sample_id}_combined"
    """
    set -euo pipefail
    total=0
    for bam in ${bam_list.join(' ')}; do
        n=\$(\$SAMTOOLS view -c "\$bam" 2>/dev/null || echo 0)
        total=\$((total + n))
    done
    if [ "\$total" -gt 0 ]; then
        \$SAMTOOLS merge -@ ${task.cpus} ${outname}.unsorted.bam ${bam_list.join(' ')}
        \$SAMTOOLS sort  -@ ${task.cpus} -o ${outname}.bam ${outname}.unsorted.bam
        \$SAMTOOLS index ${outname}.bam
        rm -f ${outname}.unsorted.bam
    else
        echo "[INFO] All input BAMs for ${sample_id} empty — creating empty combined BAM"
        printf '@HD\\tVN:1.6\\tSO:coordinate\\n' | \$SAMTOOLS view -bS -o ${outname}.bam -
        \$SAMTOOLS index ${outname}.bam
    fi
    """
}
