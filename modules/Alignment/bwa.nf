include { reference_error ; amr_error ; annotation_error } from "$baseDir/modules/nf-functions.nf"

if( params.reference ) {
    reference = file(params.reference)
    if( !reference.exists() ) return reference_error(reference)
}
if( params.amr ) {
    amr = file(params.amr)
    if( !amr.exists() ) return amr_error(amr)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}
threads = params.threads

process index {
    label "alignment"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/BWA_Index", mode: "copy"

    input:
    path fasta

    output: 
    path("${fasta}*"), emit: bwaindex

    script:
    """
    bwa index ${fasta}
    #--threads $task.cpus 
    """
}


process bwa_align {
    tag "$pair_id"
    label "alignment"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Alignment/SAM_files", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".alignment.sam") > 0) "Standard/$filename"
            else if(filename.indexOf(".alignment.dedup.sam") > 0) "Deduped/$filename"
            else {}
        }

    input:
        path dbfasta
        path indexfiles 
        tuple val(pair_id), path(reads) 

    output:
        tuple val(pair_id), path("${pair_id}.alignment.dedup.bam"), emit: bwa_dedup_bam
        tuple val(pair_id), path("${pair_id}.alignment.sorted.bam"), emit: bwa_bam

    """
     ${BWA} mem ${dbfasta} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}.alignment.sam
     ${SAMTOOLS} view -S -b ${pair_id}.alignment.sam > ${pair_id}.alignment.bam
     rm ${pair_id}.alignment.sam
     ${SAMTOOLS} sort -n ${pair_id}.alignment.bam -o ${pair_id}.alignment.sorted.bam
     rm ${pair_id}.alignment.bam
     ${SAMTOOLS} fixmate ${pair_id}.alignment.sorted.bam ${pair_id}.alignment.sorted.fix.bam
     ${SAMTOOLS} sort ${pair_id}.alignment.sorted.fix.bam -o ${pair_id}.alignment.sorted.fix.sorted.bam
     rm ${pair_id}.alignment.sorted.fix.bam 
     ${SAMTOOLS} rmdup -S ${pair_id}.alignment.sorted.fix.sorted.bam ${pair_id}.alignment.dedup.bam
     rm ${pair_id}.alignment.sorted.fix.sorted.bam
     ${SAMTOOLS} view -h -o ${pair_id}.alignment.dedup.sam ${pair_id}.alignment.dedup.bam
     rm ${pair_id}.alignment.dedup.sam 
    """
}

process bwa_rm_contaminant_fq {
    tag { pair_id }
    label "alignment"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3 
 
    publishDir "${params.output}/HostRemoval", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("fastq.gz") > 0) "NonHostFastq/$filename"
            else {}
        }

    input:
    path hostfasta
    path indexes
    tuple val(pair_id), path(reads) 

    output:
    tuple val(pair_id), path("${pair_id}.non.host.R*.fastq.gz"), emit: nonhost_reads
    path("${pair_id}.samtools.idxstats"), emit: host_rm_stats
    
    """
    ${BWA} mem ${hostfasta} ${reads[0]} ${reads[1]} -t ${threads} > ${pair_id}.host.sam
    ${SAMTOOLS} view -bS ${pair_id}.host.sam | ${SAMTOOLS} sort -@ ${threads} -o ${pair_id}.host.sorted.bam
    rm ${pair_id}.host.sam
    ${SAMTOOLS} index ${pair_id}.host.sorted.bam && ${SAMTOOLS} idxstats ${pair_id}.host.sorted.bam > ${pair_id}.samtools.idxstats
    ${SAMTOOLS} view -h -f 4 -b ${pair_id}.host.sorted.bam -o ${pair_id}.host.sorted.removed.bam
    ${BEDTOOLS}  \
       bamtofastq \
      -i ${pair_id}.host.sorted.removed.bam \
      -fq ${pair_id}.non.host.R1.fastq.gz \
      -fq2 ${pair_id}.non.host.R2.fastq.gz

    rm *.host.sam
    rm *.bam
    """

}

process HostRemovalStats {
    tag { sample_id }
    label "alignment"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
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