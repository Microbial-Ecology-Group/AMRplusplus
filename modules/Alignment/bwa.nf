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


process index {
    //tag "$referenceindex.simpleName"
    publishDir "${params.output}/BuildBWAIndex", mode: "copy"
    conda = "$baseDir/envs/AMR++_alignment.yaml"
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

    conda = "$baseDir/envs/AMR++_alignment.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'

    input:
    path dbfasta
    path indexfiles 
    tuple val(pair_id), path(reads) 

    output:
    tuple val(pair_id), path("${pair_id}.bam"), emit: bwa_bam
    tuple val(pair_id), path("${pair_id}.sam"), emit: bwa_sam

    script:
    """
    bwa mem -t ${task.cpus} ${dbfasta} ${reads} > ${pair_id}.sam
    samtools view -@ ${task.cpus} -Sb > ${pair_id}.bam
    # ${params.EXTRAPARS} 

    """
}

process bwa_rm_contaminant_fq {
    tag { pair_id }

    publishDir "${params.output}/AlignReadsToHost", mode: "copy"

    conda = "$baseDir/envs/AMR++_alignment.yaml"
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