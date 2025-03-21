// Resistome

if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}

threshold = params.threshold
threads = params.threads


min = params.min
max = params.max
skip = params.skip
samples = params.samples

deduped = params.deduped
prefix = params.prefix

process build_dependencies {
    tag { dl_dependencies }
    label "python"

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
    samtools view -h ${bam} > ${sample_id}.sam
    
    $resistome -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sample_id}.sam \
      -gene_fp ${sample_id}.${prefix}.gene.tsv \
      -group_fp ${sample_id}.${prefix}.group.tsv \
      -mech_fp ${sample_id}.${prefix}.mechanism.tsv \
      -class_fp ${sample_id}.${prefix}.class.tsv \
      -type_fp ${sample_id}.${prefix}.type.tsv \
      -t ${threshold}

    rm ${sample_id}.sam
    """
}

process resistomeresults {
    tag "Make AMR count matrix"
    label "python"

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
    samtools view -h ${bam} > ${sample_id}.sam

    $rarefaction \
      -ref_fp ${amr} \
      -sam_fp ${sample_id}.sam \
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

    rm ${sample_id}.sam
    """
}

process plotrarefaction {
    tag "Plot rarefaction results"
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".png") > 0) "Rarefaction/Figures/$filename"
            else {}
        }

    input:
        path(rarefaction)

    output:
        path("*.png"), emit: plots

    """
    mkdir -p data/
    mv *.tsv data/
    python $baseDir/bin/rfplot.py --dir ./data --nd --s --sd .
    """
}


process runsnp {
    tag {sample_id}
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/ResistomeAnalysis/SNP_verification", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".csv") > 0) "SNP_raw_output/$filename"
            else if(filename.indexOf(".txt") > 0 ) "SNP_resistant_reads/$filename"
            else if(filename.indexOf(".tsv") > 0 ) "SNP_verification_counts/$filename"
            else {}
        }

    errorStrategy = 'ignore'

    input:
        tuple val(sample_id), path(bam)
        path(snp_count_matrix)

    output:
        path("${sample_id}.SNP_confirmed_gene.tsv"), emit: snp_counts
        path("${sample_id}.${prefix}_SNPs${sample_id}/*")
        path("${sample_id}_${prefix}_SNPresistant_reads.txt")

    """
    cp -r $baseDir/bin/AmrPlusPlus_SNP/* .

    # change name to stay consistent with count matrix name, but only if the names don't match
    if [ "${bam}" != "${sample_id}.bam" ]; then
        mv ${bam} ${sample_id}.bam
    fi

    python3 SNP_Verification.py -c config.ini -t ${threads} -a true -i ${sample_id}.bam -o ${sample_id}.${prefix}_SNPs --count_matrix ${snp_count_matrix} --detailed_output = True

    cut -d ',' -f `awk -v RS=',' "/${sample_id}/{print NR; exit}" ${sample_id}.${prefix}_SNPs${snp_count_matrix}` ${sample_id}.${prefix}_SNPs${snp_count_matrix} > ${sample_id}.${prefix}_SNP_count_col

    cut -d ',' -f 1 ${sample_id}.${prefix}_SNPs${snp_count_matrix} > gene_accession_labels

    paste gene_accession_labels ${sample_id}.${prefix}_SNP_count_col > ${sample_id}.SNP_confirmed_gene.tsv

    mv ${sample_id}.${prefix}_SNPs${sample_id}/resistant_reads.csv ${sample_id}_${prefix}_SNPresistant_reads.txt

    """
}


process snpresults {
    tag "Make SNP-confirmed matrix"
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    publishDir "${params.output}/Results", mode: "copy"

    errorStrategy = 'ignore'

    input:
        path(snp_counts)

    output:
        path("*_analytic_matrix.csv"), emit: snp_matrix

    """

    ${PYTHON3} $baseDir/bin/snp_long_to_wide.py -i ${snp_counts} -o SNPconfirmed_${prefix}_analytic_matrix.csv

    """
}
