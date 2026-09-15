// modules/SNV/metasnv.nf
//
// SNV calling from BAM files aligned to MEGARes, following the
// NGLess -> metaSNV -> parsing workflow.
//
// ═══════════════════════════════════════════════════════════════════════════
// THE REFERENCE MUST MATCH THE BAMs
// ═══════════════════════════════════════════════════════════════════════════
// metaSNV associates reads with reference sequences by contig name, so the BAMs
// supplied here must have been aligned to the same fasta passed as
// --snv_reference. A mismatch fails immediately or, worse, produces empty
// results with no obvious cause.
//
// That reference should normally be the MEGARes REPRESENTATIVE database. Both
// the NGLess filter below and metaSNV itself depend on reads mapping to exactly
// one reference sequence, because an allele cannot be attributed to a gene when
// the read maps equally well to several. The complete database contains large
// families of near-identical accessions, so a large share of short reads become
// multi-mappers and are discarded by the unique-mapping filter.
// ═══════════════════════════════════════════════════════════════════════════


// ─────────────────────────────────────────────────────────────────────────────
// NGLess filtering, one sample at a time.
//
// The .ngl script is generated here rather than by bin/SNV_calling/prep_NGLess.py
// for two reasons. First, prep_NGLess.py emits a script built around
// run_for_all(readlines(...)), which makes NGLess iterate over every sample
// itself; that duplicates the parallelism Nextflow already provides and prevents
// per-sample resume and retry. Second, it writes the script to bin/filter.ngl,
// which is outside the task work directory and read-only under containers.
//
// prep_NGLess.py remains useful for running the filter standalone outside the
// pipeline. The filter criteria below are identical to the ones it generates.
// ─────────────────────────────────────────────────────────────────────────────
process ngless_filter_snv {
    tag { sample_id }
    label "medium"

    publishDir "${params.output}/SNV_analysis/FilteredBAMs", mode: "copy"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.SNV.unique.sorted.bam"), emit: filtered_bam

    script:
    """
    set -euo pipefail

    cat > filter_${sample_id}.ngl <<'NGLEOF'
        ngless "1.5"
        import "samtools" version "0.0"

        # Filter alignments for SNV calling with metaSNV v2.
        # Keep only alignments meeting the minimum aligned-block size and percent
        # identity, then keep only reads that mapped uniquely. Multi-mapped reads are
        # discarded because a variant cannot be attributed to one gene when the read
        # maps equally well to another.

        input = samfile('${bam}')

        filtered = select(input) using |mr|:
            mr = mr.filter(min_match_size=${params.snv_min_match_size}, min_identity_pc=${params.snv_min_identity_pc}, action={unmatch})

        filtered_unique = select(filtered, keep_if=[{mapped}, {unique}])
        filtered_unique = samtools_sort(filtered_unique)
        write(filtered_unique, ofile='${sample_id}.SNV.unique.sorted.bam')
        NGLEOF

    ngless --trace -j ${task.cpus} filter_${sample_id}.ngl
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Index the filtered BAMs.
//
// NGLess samtools_sort produces coordinate-sorted output, so no re-sort is
// needed, but metaSNV requires an index alongside each BAM.
// ─────────────────────────────────────────────────────────────────────────────
process index_snv_bams {
    tag { sample_id }
    label "small"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    set -euo pipefail
    \$SAMTOOLS index -@ ${task.cpus} ${bam}
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// metaSNV, run ONCE across all samples.
//
// metaSNV derives per-position statistics across the whole sample set, so it is
// not a per-sample process. Two staging details matter: metaSNV writes index and
// gen_pos files next to the reference, so the fasta is copied rather than used
// through the staged symlink; and 'all_samples' must list one BAM path per line,
// in the order that becomes the column order of the count matrix.
// ─────────────────────────────────────────────────────────────────────────────
process run_metasnv {
    tag "metaSNV"
    label "large"

    publishDir "${params.output}", mode: "copy"

    input:
        path(bams)
        path(bais)
        path(reference)

    output:
        path("SNV_analysis_output"),             emit: snv_dir
        path("SNV_analysis_output/all_samples"), emit: all_samples

    script:
    def db_ann   = params.snv_db_ann   ? "--db_ann ${params.snv_db_ann}"     : ""
    def n_splits = params.snv_n_splits ? "--n_splits ${params.snv_n_splits}" : ""
    """
    set -euo pipefail

    # metaSNV writes auxiliary files beside the reference, so give it a real
    # writable copy rather than the staged symlink.
    cp -L ${reference} snv_reference.fasta
    chmod u+w snv_reference.fasta

    # One BAM path per line. Sorted so the column order of the resulting matrix
    # is deterministic across runs rather than depending on channel arrival order.
    for b in \$(ls -1 *.bam | sort); do
        n=\$(\$SAMTOOLS view -c "\$b" 2>/dev/null || echo 0)
        if [ "\$n" -gt 0 ]; then
            readlink -f "\$b" >> all_samples
        else
            echo "[WARN] skipping \$b: no alignments survived filtering" >&2
        fi
    done

    if [ ! -s all_samples ]; then
        echo "[ERROR] No non-empty BAM files remain after NGLess filtering." >&2
        echo "[ERROR] Check snv_min_match_size and snv_min_identity_pc." >&2
        exit 1
    fi

    echo "[INFO] Running metaSNV on \$(wc -l < all_samples) sample(s)"

    metaSNV.py SNV_analysis_output all_samples snv_reference.fasta \\
        --threads ${task.cpus} ${n_splits} ${db_ann}

    # clean_metaSNP reads the run-order sample list from the output directory.
    if [ ! -f SNV_analysis_output/all_samples ]; then
        cp all_samples SNV_analysis_output/all_samples
    fi
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Reformat raw metaSNV calls into an AMR++-style count matrix and annotation
// file. Wraps bin/SNV_calling/clean_metaSNP_1.0.2.py.
//
// That script drops RequiresSNPConfirmation genes, which are handled by the
// separate SNP confirmation workflow, and builds one SNV_accession per variant:
//     MEG_id|Type|Class|Mechanism|Group|GROUP-accession_REF_POS_ALT
// so individual variants can be analysed as distinct resistome features.
// ─────────────────────────────────────────────────────────────────────────────
process clean_metasnv {
    tag "clean_metaSNV"
    label "medium"

    publishDir "${params.output}/Results", mode: "copy"

    input:
        path(snv_dir)

    output:
        path("*SNV_analytic_matrix.csv"), emit: snv_matrix
        path("*SNV*annotations.csv"),     emit: snv_annotations

    script:
    def wf = (params.snv_aln_wf == "Deduped") ? "Deduped" : "Standard"
    def custom_matrix = params.snv_matrix_out ? "-matrix_out ${params.snv_matrix_out}" : ""
    def custom_ann    = params.snv_ann_out    ? "-ann_out ${params.snv_ann_out}"       : ""
    """
    set -euo pipefail

    RAW_DIR="${snv_dir}/snpCaller"
    if [ ! -d "\$RAW_DIR" ]; then
        echo "[ERROR] \$RAW_DIR not found; metaSNV produced no raw SNP calls." >&2
        ls -la ${snv_dir} >&2 || true
        exit 1
    fi

    \$PYTHON3 $baseDir/bin/SNV_calling/clean_metaSNP_1.0.2.py \\
        -s ./ \\
        -file_dir "\$RAW_DIR/" \\
        -fn_struc "${params.snv_snpfile_prefix}" \\
        -sample_file "${snv_dir}/all_samples" \\
        -aln_wf ${wf} \\
        ${custom_matrix} \\
        ${custom_ann}
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// Restrict the SNV matrix to genes also detected in the resistome count matrix.
// Wraps bin/SNV_calling/Filter_SNV_matrix.py.
//
// Requires the resistome matrix, so this only runs when the resistome workflow
// produced one. Optional; controlled by params.snv_filter_by_resistome.
// ─────────────────────────────────────────────────────────────────────────────
process filter_snv_matrix {
    tag "filter_SNV"
    label "medium"

    publishDir "${params.output}/Results", mode: "copy"

    input:
        path(snv_matrix)
        path(resistome_matrix)

    output:
        path("filtered_*.csv"), emit: filtered_matrix

    script:
    """
    set -euo pipefail

    \$PYTHON3 $baseDir/bin/SNV_calling/Filter_SNV_matrix.py \\
        --snv-matrix ${snv_matrix} \\
        --resistome-matrix ${resistome_matrix} \\
        --output filtered_${snv_matrix}
    """
}
