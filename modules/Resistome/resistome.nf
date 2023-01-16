// Resistome

if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}

threshold = params.threshold

min = params.min
max = params.max
skip = params.skip
samples = params.samples

deduped = params.deduped
prefix = params.prefix

process build_dependencies {
    tag { dl_dependencies }
    label "python"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${baseDir}/bin/", mode: "copy"

    output:
        path("rarefaction"), emit: rarefactionanalyzer
        path("resistome"), emit: resistomeanalyzer
        path("AmrPlusPlus_SNP/*"), emit: amrsnp

    """
    # Uncomment these sections once the v2 rarefactionanalyzer and resistomeanalyzer repositories are updated, remove cp lines
    #git clone https://github.com/cdeanj/rarefactionanalyzer.git
    #cd rarefactionanalyzer
    #make
    #chmod 777 rarefaction
    #mv rarefaction ../
    #cd ../
    #rm -rf rarefactionanalyzer
    cp $baseDir/bin/rarefaction . 


    #git clone https://github.com/cdeanj/resistomeanalyzer.git
    #cd resistomeanalyzer
    #make
    #chmod 777 resistome
    #mv resistome ../
    #cd ../
    #rm -rf resistomeanalyzer
    cp $baseDir/bin/resistome .

    git clone https://github.com/Isabella136/AmrPlusPlus_SNP.git
    chmod -R 777 AmrPlusPlus_SNP/

    """


}


process runresistome {
    tag { sample_id }
    label "alignment"

    memory { 3.GB * task.attempt }
    time { 6.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".tsv") > 0) "ResistomeCounts/$filename"
            else {}
        }

    input:
        tuple val(sample_id), path(bam)
        path(amr)
        path(annotation)
        path(resistome)

    output:
        tuple val(sample_id), path("${sample_id}*.tsv"), emit: resistome_tsv
        path("${sample_id}.${prefix}.gene.tsv"), emit: resistome_counts

    
    
    """
    samtools view -h ${bam} > converted.sam
    
    $resistome -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp converted.sam \
      -gene_fp ${sample_id}.${prefix}.gene.tsv \
      -group_fp ${sample_id}.${prefix}.group.tsv \
      -mech_fp ${sample_id}.${prefix}.mechanism.tsv \
      -class_fp ${sample_id}.${prefix}.class.tsv \
      -type_fp ${sample_id}.${prefix}.type.tsv \
      -t ${threshold}

    rm converted.sam
    """
}

process resistomeresults {
    tag { }
    label "python"

    memory { 3.GB * task.attempt }
    time { 4.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/Results", mode: "copy"

    input:
        path(resistomes)

    output:
        path("${prefix}_analytic_matrix.csv"), emit: raw_count_matrix
        path("${prefix}_analytic_matrix.csv"), emit: snp_count_matrix, optional: true

    """
    ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o ${prefix}_analytic_matrix.csv
    """
}

process runrarefaction {
    tag { sample_id }
    label "alignment"

    memory { 3.GB * task.attempt }
    time { 3.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".tsv") > 0) "Rarefaction/Counts/$filename"
            else {}
        }

    input:
        tuple val(sample_id), path(bam)
        path(annotation)
        path(amr)
        path(rarefaction)

    output:
        path("*.tsv"), emit: rarefaction

    """
    samtools view -h ${bam} > converted.sam

    $rarefaction \
      -ref_fp ${amr} \
      -sam_fp converted.sam \
      -annot_fp ${annotation} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -mech_fp ${sample_id}.mech.tsv \
      -class_fp ${sample_id}.class.tsv \
      -type_fp ${sample_id}.type.tsv \
      -min ${min} \
      -max ${max} \
      -skip ${skip} \
      -samples ${samples} \
      -t ${threshold}

    rm converted.sam
    """
}

process plotrarefaction {
    tag { sample_id }
    label "python"

    memory { 3.GB * task.attempt }
    time { 2.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("graphs/*.png") > 0) "Rarefaction/Figures/$filename"
            else {}
        }

    input:
        path(rarefaction)

    output:
        path("graphs/*.png"), emit: plots

    """
    mkdir data/
    mv *.tsv data/
    mkdir graphs/
    python $baseDir/bin/rfplot.py --dir ./data --nd --s --sd ./graphs
    """
}


process runsnp {
    tag {sample_id}
    label "python"


    memory { 3.GB * task.attempt }
    time { 2.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf("_SNPs/*") > 0) "SNP_verification/$filename"
            else {}
        }

    errorStrategy = 'ignore'

    input:
        tuple val(sample_id), path(bam)
        path(snp_count_matrix)

    output:
        path("${sample_id}*_count_col"), emit: snp_counts
        path("${sample_id}*_SNPs/*")

    """
    cp -r $baseDir/bin/AmrPlusPlus_SNP/* .

    samtools view -h ${bam} > converted.sam

    python3 SNP_Verification.py -c config.ini -a -i converted.sam -o ${sample_id}_${prefix}_SNPs --count_matrix ${snp_count_matrix}

    cut -d ',' -f `awk -v RS=',' "/${sample_id}/{print NR; exit}" ${snp_count_matrix}` ${snp_count_matrix} > ${sample_id}_${prefix}_SNP_count_col

    rm converted.sam

    """
}


process snpresults {
    tag {sample_id}
    label "python"

    memory { 3.GB * task.attempt }
    time { 2.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/Results", mode: "copy"

    errorStrategy = 'ignore'

    input:
        path(snp_counts)
        path(snp_count_matrix)

    output:
        path("*_analytic_matrix.csv"), emit: snp_matrix

    """

    cut -d ',' -f 1 ${snp_count_matrix} > gene_accession_labels
    paste gene_accession_labels ${snp_counts} > SNPconfirmed_${prefix}_analytic_matrix.csv


    """
}
