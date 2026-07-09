// Resistome

process build_dependencies {
    tag "Download SNP dependencies"
    label "nano"

    publishDir "${baseDir}/bin/", mode: "copy"

    output:
        path("rarefaction"), emit: rarefactionanalyzer
        path("resistome"), emit: resistomeanalyzer
        path("AmrPlusPlus_SNP/*"), emit: amrsnp

    script:
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

    #git clone https://github.com/Isabella136/AmrPlusPlus_SNP.git
    git clone https://github.com/EnriqueDoster/AmrPlusPlus_SNP.git
    chmod -R 777 AmrPlusPlus_SNP/
    cd AmrPlusPlus_SNP
    """


}


process runresistome {
    tag { sample_id }
    label "medium"

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
        path("${sample_id}.${params.prefix}.gene.tsv"), emit: resistome_counts

    script:
    """
    set -euo pipefail

    # Check if BAM has any alignments
    aln_count=\$(\$SAMTOOLS view -c ${bam} 2>/dev/null || echo 0)

    if [ "\$aln_count" -gt 0 ]; then
        \$SAMTOOLS view -h ${bam} > ${sample_id}.sam

        $resistome -ref_fp ${amr} \\
          -annot_fp ${annotation} \\
          -sam_fp ${sample_id}.sam \\
          -gene_fp ${sample_id}.${params.prefix}.gene.tsv \\
          -group_fp ${sample_id}.${params.prefix}.group.tsv \\
          -mech_fp ${sample_id}.${params.prefix}.mechanism.tsv \\
          -class_fp ${sample_id}.${params.prefix}.class.tsv \\
          -type_fp ${sample_id}.${params.prefix}.type.tsv \\
          -t ${params.threshold}

        rm ${sample_id}.sam
    else
        echo "[INFO] No alignments in BAM for ${sample_id} — writing empty resistome counts"
        for level in gene group mechanism class type; do
            printf "Header\\t0\\n" > ${sample_id}.${params.prefix}.\${level}.tsv
        done
    fi
    """
}


process runresistome_analyzer {
    tag { sample_id }
    label "medium"

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if      (filename.endsWith(".gene.tsv"))            "ResistomeCounts/$filename"
            else if (filename.endsWith("_per_read.tsv"))        "PerReadInfo/$filename"
            else if (filename.endsWith("_coverage_stats.tsv"))  "CoverageStats/$filename"
            else {}
        }

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}*.tsv"),    emit: resistome_tsv
        path("${sample_id}.${params.prefix}.gene.tsv"),      emit: resistome_counts

    script:
    def count_mode = params.count_mode ?: "fragment"
    def group_flag = (params.group_aware == "N")     ? "--no-group-aware"     : "--group-aware"
    def edge_flag  = (params.edge_aware_qcov == "N")  ? "--no-edge-aware-qcov" : "--edge-aware-qcov"
    def supp_flag  = (params.include_supplementary == "Y") ? "--include-supplementary" : ""
    def cigar_flag = (params.cigar_aware_coverage  == "Y") ? "--cigar-aware-coverage"  : ""
    """
    set -euo pipefail

    aln_count=\$(\$SAMTOOLS view -c ${bam} 2>/dev/null || echo 0)

    if [ "\$aln_count" -gt 0 ]; then
        \$PYTHON3 $baseDir/bin/alignment_analyzer.py \\
            -i ${bam} \\
            -r ${sample_id}.${params.prefix}_per_read.tsv \\
            -g ${sample_id}.${params.prefix}.gene_summary.tsv \\
            --sample-id ${sample_id} \\
            --count-mode ${count_mode} \\
            --min-mapq ${params.min_mapq} \\
            --min-gene-fraction ${params.min_gene_fraction} \\
            --min-query-coverage ${params.min_query_coverage} \\
            --coverage-output ${sample_id}.${params.prefix}_coverage_stats.tsv \\
            ${group_flag} ${edge_flag} ${supp_flag} ${cigar_flag}

        # Reshape: gene_accession \\t meg_id \\t count  ->  sample \\t gene \\t count
        tail -n +2 ${sample_id}.${params.prefix}.gene_summary.tsv \\
          | awk -v s="${sample_id}" 'BEGIN{FS=OFS="\\t"} {print s, \$1, \$3}' \\
          > ${sample_id}.${params.prefix}.gene.tsv
    else
        echo "[INFO] No alignments in BAM for ${sample_id} — writing empty gene counts"
        : > ${sample_id}.${params.prefix}.gene.tsv
    fi
    """
}

// Coverage-threshold sweep diagnostic. Runs coverage_threshold_sweep.py on each BAM
// to profile how gene detection drops off as coverage/quality thresholds increase,
// then summarizes all per-sample sweeps with plot_sweep_dropoff.R.
//
// NOTE on counting: coverage_threshold_sweep.py counts per-FRAGMENT (mates resolved
// together via base_read_id), so paired-end (unmerged) and FLASH-merged reads are
// weighted identically. This matches --count-mode fragment in alignment_analyzer.py.

process coverage_threshold_sweep {
    tag { prefix }
    label "medium"
    publishDir "${params.output}/CoverageSweep/PerSample", mode: "copy"

    input:
        tuple val(prefix), path(bam)

    output:
        // All four are emitted into one channel so the combine step can stage them
        // together. gene_detail + length_quantiles are what plot_sweep_dropoff.R needs;
        // _results.csv is the main grid; _redundancy.csv is optional diagnostics.
        path("${prefix}_*.csv"), emit: sweep_csvs

    script:
    def group_flag = (params.group_aware == "N")          ? "--no-group-aware"     : "--group-aware"
    def edge_flag  = (params.sweep_edge_aware_qcov == "N") ? "--no-edge-aware-qcov" : "--edge-aware-qcov"
    def snp_flag   = (params.sweep_exclude_snp == "Y")     ? "--exclude-snp-confirmation" : ""
    """
    set -euo pipefail

    aln_count=\$(\$SAMTOOLS view -c ${bam} 2>/dev/null || echo 0)

    if [ "\$aln_count" -gt 0 ]; then
        \$PYTHON3 $baseDir/bin/coverage_threshold_sweep.py \\
            -i ${bam} \\
            -o ${prefix}_results.csv \\
            --min-mapq ${params.min_mapq} \\
            --gene-detail-output ${prefix}_gene_detail.csv \\
            --redundancy-output ${prefix}_redundancy.csv \\
            --length-quantiles-output ${prefix}_length_quantiles.csv \\
            ${group_flag} ${edge_flag} ${snp_flag}
    else
        echo "[INFO] No alignments in ${bam} — writing empty sweep CSVs"
        echo "min_query_coverage,min_identity,min_match_qcov,min_gene_fraction,n_alignments_retained" \\
            > ${prefix}_results.csv
        echo "min_query_coverage,gene_accession,gene_fraction,read_count" > ${prefix}_gene_detail.csv
        echo "min_query_coverage,group,n_accessions" > ${prefix}_redundancy.csv
        echo "min_query_coverage,n_alignments" > ${prefix}_length_quantiles.csv
    fi
    """
}

process plot_sweep_dropoff {
    tag "sweep_summary"
    label "small"
    publishDir "${params.output}/CoverageSweep/Summary", mode: "copy"

    input:
        path(sweep_csvs)

    output:
        path("*.pdf"), optional: true
        path("*.png"), optional: true
        path("*.tsv"), optional: true
        path("*.csv"), optional: true

    script:
    """
    set -euo pipefail
    \$RSCRIPT $baseDir/bin/plot_sweep_dropoff.R ${sweep_csvs.join(' ')}
    """
}


process resistomeresults {
    tag "Make AMR count matrix"
    label "small"
    
    publishDir "${params.output}/Results", mode: "copy"

    input:
        path(resistomes)
        val  prefix

    output:
        path("${params.prefix}_analytic_matrix.csv"), emit: raw_count_matrix
        path("${params.prefix}_analytic_matrix.csv"), emit: snp_count_matrix, optional: true
    script:
    """
    \$PYTHON3 $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o ${params.prefix}_analytic_matrix.csv
    """
}

process runrarefaction {
    tag { sample_id }
    label "small"

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

    script:
    """
    set -euo pipefail

    aln_count=\$(\$SAMTOOLS view -c ${bam} 2>/dev/null || echo 0)

    if [ "\$aln_count" -gt 0 ]; then
        \$SAMTOOLS view -h ${bam} > ${sample_id}.sam

        $rarefaction \\
          -ref_fp ${amr} \\
          -sam_fp ${sample_id}.sam \\
          -annot_fp ${annotation} \\
          -gene_fp ${sample_id}.gene.tsv \\
          -group_fp ${sample_id}.group.tsv \\
          -mech_fp ${sample_id}.mech.tsv \\
          -class_fp ${sample_id}.class.tsv \\
          -type_fp ${sample_id}.type.tsv \\
          -min ${params.min} \\
          -max ${params.max} \\
          -skip ${params.skip} \\
          -samples ${params.samples} \\
          -t ${params.threshold}

        rm ${sample_id}.sam
    else
        echo "[INFO] No alignments in BAM for ${sample_id} — writing empty rarefaction counts"
        for level in gene group mech class type; do
            printf "0\\t0\\n" > ${sample_id}.\${level}.tsv
        done
    fi
    """
}


process plotrarefaction {
    tag "Plot rarefaction results"
    label "micro"

    publishDir "${params.output}/ResistomeAnalysis", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".png") > 0) "Rarefaction/Figures/$filename"
            else {}
        }

    input:
        path(rarefaction)
        val  prefix

    output:
        path("*.png"), emit: plots

    script:
    """
    mkdir -p data/
    mv *.tsv data/
    python $baseDir/bin/rfplot.py --dir ./data --nd --s --sd . --prefix ${prefix}
    """
}


process old_runsnp {
    tag {sample_id}
    label "snp_ignore"

    publishDir "${params.output}/ResistomeAnalysis/SNP_verification", mode: "copy",
            saveAs: { filename ->
                if(filename.indexOf(".tsv") > 0) "SNP_verification_counts/$filename"
                else "SNP_detailed_output/$filename"
            }

    input:
        tuple val(sample_id), path(bam)
        path(snp_count_matrix)

    output:
        path("${sample_id}.SNP_confirmed_gene.tsv"), emit: snp_counts
    script:
    """
    cp -rsa $baseDir/bin/AmrPlusPlus_SNP/* .

    # change name to stay consistent with count matrix name, but only if the names don't match
    if [ "${bam}" != "${sample_id}.bam" ]; then
        mv ${bam} ${sample_id}.bam
    fi

    \$PYTHON3 SNP_Verification.py -c config.ini -t ${task.cpus} -a true -i ${sample_id}.bam -o ${sample_id}.${params.prefix}_SNPs --count_matrix ${snp_count_matrix} --detailed_output=all

    \$PYTHON3 $baseDir/bin/extract_snp_column.py \
      --sample-id "${sample_id}" \
      --matrix "${sample_id}.${params.prefix}_SNPs${snp_count_matrix}" \
      --out-tsv "${sample_id}.SNP_confirmed_gene.tsv"
    """
}

process runsnp {
    tag {sample_id}
    label "snp_ignore"

    publishDir "${params.output}/ResistomeAnalysis/SNP_verification", mode: "copy",
            saveAs: { filename ->
                if(filename.indexOf(".tsv") > 0) "SNP_verification_counts/$filename"
                else "SNP_detailed_output/$filename"
            }

    input:
        tuple val(sample_id), path(bam)
        path(snp_count_matrix)

    output:
        path("${sample_id}.SNP_confirmed_gene.tsv"), emit: snp_counts
        path("${sample_id}.${params.prefix}_SNPs/${sample_id}/${sample_id}.${params.prefix}_SNPs_SNPs_resistant_reads.txt") , optional: true
        path("${sample_id}.${params.prefix}_SNPs/${sample_id}/${sample_id}.${params.prefix}_SNPs_snp_coverage_stats.csv")
        path("${sample_id}.${params.prefix}_SNPs/${sample_id}/${sample_id}.${params.prefix}_SNPs_snp_verification_summary.csv")
    script:
    """
    cp -rsa $baseDir/bin/AmrPlusPlus_SNP/* .

    # change name to stay consistent with count matrix name, but only if the names don't match
    if [ "${bam}" != "${sample_id}.bam" ]; then
        mv ${bam} ${sample_id}.bam
    fi

    \$PYTHON3 SNP_Verification.py -c config.ini -t ${task.cpus} -a true -i ${sample_id}.bam -o ${sample_id}.${params.prefix}_SNPs --count_matrix ${snp_count_matrix} --detailed_output=all

    \$PYTHON3 $baseDir/bin/extract_snp_column.py \
      --sample-id "${sample_id}" \
      --matrix ${sample_id}.${params.prefix}_SNPs/"${sample_id}.${params.prefix}_SNPs_${snp_count_matrix}" \
      --out-tsv "${sample_id}.SNP_confirmed_gene.tsv"
    """
}

process snpresults {
    tag "Make SNP-confirmed matrix"
    label "micro"
    
    publishDir "${params.output}/Results", mode: "copy"

    input:
        path(snp_counts)
        val  prefix

    output:
        path("*_analytic_matrix.csv"), emit: snp_matrix

    script:
    """
    \$PYTHON3 $baseDir/bin/snp_long_to_wide.py -i ${snp_counts} -o SNPconfirmed_${prefix}_analytic_matrix.csv
    """
}


