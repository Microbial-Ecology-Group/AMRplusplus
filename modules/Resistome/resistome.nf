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

process runresistome {
    tag { sample_id }

    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        tuple val(sample_id), path(sam)
        path(amr)
        path(annotation)

    output:
        tuple val(sample_id), path("${sample_id}*.tsv"), emit: resistome_tsv
        path("${sample_id}.gene.tsv"), emit: resistome_counts

    """
    $baseDir/bin/resistome -ref_fp ${amr} \
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

    publishDir "${params.output}/RunSNP_Verification", mode: "copy"

    input:
        tuple val(sample_id), path(sam)

    output:
        tuple val(sample_id), path("${sample_id}_SNPs/*"), emit: snps

    """
    python3 $baseDir/AmrPlusPlus_SNP/SNP_Verification.py -i ${sam} -o ${sample_id}_SNPs
    """
}

process resistomeresults {
    tag { }

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

    publishDir "${params.output}/RunRarefaction", mode: "copy"

    input:
        tuple val(sample_id), path(sam)
        path(annotation)
        path(amr)

    output:
        path("*.tsv"), emit: rarefaction

    """
    $baseDir/bin/rarefaction \
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