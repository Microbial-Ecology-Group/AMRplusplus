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
    //tag "$referenceindex.simpleName"
    publishDir "${params.output}/BuildBWAIndex", mode: "copy"
    conda = "$baseDir/envs/alignment.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'

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
    publishDir "${params.output}/AlignToDB", mode: "copy"

    conda = "$baseDir/envs/alignment.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'

    input:
    path dbfasta
    path indexfiles 
    tuple val(pair_id), path(reads) 

    output:
    tuple val(pair_id), path("${pair_id}.amr.alignment.dedup.bam"), emit: bwa_dedup_bam
    tuple val(pair_id), path("${pair_id}.amr.alignment.sorted.fix.sorted.bam"), emit: bwa_bam
    tuple val(pair_id), path("${pair_id}.amr.alignment.dedup.sam"), emit: bwa_dedup_sam
    tuple val(pair_id), path("${pair_id}.amr.alignment.sam"), emit: bwa_sam

    script:
    """
     ${BWA} mem ${dbfasta} ${reads} -t ${threads} -R '@RG\\tID:${pair_id}\\tSM:${pair_id}' > ${pair_id}.amr.alignment.sam
     ${SAMTOOLS} view -S -b ${pair_id}.amr.alignment.sam > ${pair_id}.amr.alignment.bam
     ${SAMTOOLS} sort -n ${pair_id}.amr.alignment.bam -o ${pair_id}.amr.alignment.sorted.bam
     ${SAMTOOLS} fixmate ${pair_id}.amr.alignment.sorted.bam ${pair_id}.amr.alignment.sorted.fix.bam
     ${SAMTOOLS} sort ${pair_id}.amr.alignment.sorted.fix.bam -o ${pair_id}.amr.alignment.sorted.fix.sorted.bam
     ${SAMTOOLS} rmdup -S ${pair_id}.amr.alignment.sorted.fix.sorted.bam ${pair_id}.amr.alignment.dedup.bam
     ${SAMTOOLS} view -h -o ${pair_id}.amr.alignment.dedup.sam ${pair_id}.amr.alignment.dedup.bam
     #rm ${pair_id}.amr.alignment.bam
     #rm ${pair_id}.amr.alignment.sorted*.bam
    """
}

process bwa_rm_contaminant_fq {
    tag { pair_id }

    publishDir "${params.output}/AlignReadsToHost", mode: "copy"

    conda = "$baseDir/envs/alignment.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'

    input:
    path hostfasta
    path indexes
    tuple val(pair_id), path(reads) 

    output:
    tuple val(pair_id), path("${pair_id}.non.host.R*.fastq.gz"), emit: nonhost_reads
    
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

    """


}