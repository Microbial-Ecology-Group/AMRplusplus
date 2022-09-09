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

process build_dependencies {
    tag { dl_github }
    publishDir "${baseDir}/bin/", mode: "copy"
    conda = "$baseDir/envs/python.yaml"

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

    """


}






process runresistome {
    tag { sample_id }
    conda = "$baseDir/envs/alignment.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'
    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        tuple val(sample_id), path(sam)
        path(amr)
        path(annotation)
        path(resistome)

    output:
        tuple val(sample_id), path("${sample_id}*.tsv"), emit: resistome_tsv
        path("${sample_id}.gene.tsv"), emit: resistome_counts

    
    
    """
    $resistome -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sam} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -mech_fp ${sample_id}.mechanism.tsv \
      -class_fp ${sample_id}.class.tsv \
      -type_fp ${sample_id}.type.tsv \
      -t ${threshold}
    """
}

process runsnp {
    tag { sample_id }
    conda = "$baseDir/envs/python.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'
    publishDir "${params.output}/RunSNP_Verification", mode: "copy"

    errorStrategy = 'ignore'

    input:
        tuple val(sample_id), path(sam)
        path(amrsnp)

    output:
        tuple val(sample_id), path("${sample_id}_SNPs/*"), emit: snps

    """
    python3 $amrsnp/SNP_Verification.py -i ${sam} -o ${sample_id}_SNPs
    """
}

process resistomeresults {
    tag { }
    conda = "$baseDir/envs/python.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'
    publishDir "${params.output}/ResistomeResults", mode: "copy"

    input:
        path(resistomes)

    output:
        path("AMR_analytic_matrix.csv"), emit: raw_count_matrix

    """
    ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o AMR_analytic_matrix.csv
    """
}

process runrarefaction {
    tag { sample_id }
    conda = "$baseDir/envs/python.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'
    publishDir "${params.output}/RunRarefaction", mode: "copy"

    input:
        tuple val(sample_id), path(sam)
        path(annotation)
        path(amr)
        path(rarefaction)

    output:
        path("*.tsv"), emit: rarefaction

    """
    $rarefaction \
      -ref_fp ${amr} \
      -sam_fp ${sam} \
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
    """
}

process plotrarefaction {
    tag { sample_id }
    conda = "$baseDir/envs/python.yaml"
    container = 'enriquedoster/amrplusplus_alignment:latest'
    publishDir "${params.output}/RarefactionFigures", mode: "copy"

    input:
        path(rarefaction)

    output:
        path("*.tsv"), emit: rarefaction

    """
    mkdir data/
    mv *.tsv data/
    mkdir graphs/
    python $baseDir/bin/rfplot.py --dir ./data --nd --s --sd ./graphs
    """
}
