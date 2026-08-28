# Coverage Threshold Sweep: Step-by-Step Tutorial

This tutorial covers the `bam_coverage_sweep` pipeline, which takes a set of pre-aligned BAM files and evaluates the resistome across a grid of alignment filter thresholds. Instead of committing to a single cutoff, the sweep shows you how the number of detected genes, groups, mechanisms and classes drops off as you tighten each threshold, so you can choose sensible cutoffs and see how sensitive each sample is to them. It produces combined result matrices and a set of drop-off plots.

> **Before you start:** This pipeline runs on BAM files you have already produced with AMR++ (for example, the `*_alignment_sorted.bam` files under `Alignment/BAM_files/Standard/`). It also requires a few extra tools that are **not** bundled in the default AMR++ conda environment. See [Additional tools required](#additional-tools-required) below.

## Table of Contents

- [Additional tools required](#additional-tools-required)
- [What the sweep does](#what-the-sweep-does)
- [Running the sweep](#running-the-sweep)
- [The three filters](#the-three-filters)
- [Choosing what to sweep](#choosing-what-to-sweep)
- [Outputs](#outputs)
- [Interpreting the plots](#interpreting-the-plots)
- [Applying the selected filter thresholds](#applying-the-selected-filter-thresholds)

---

## Additional tools required

The sweep uses a Python analysis step and an R plotting step whose dependencies are **not included in the AMR++ conda environment**, because adding them, especially the R plotting stack, would substantially increase an already large environment. Install them separately before running the sweep.

**Python** (for `coverage_threshold_sweep.py` and `combine_sweep_results.py`):

- `pysam` for reading BAM alignments
- `duckdb` for fast aggregation of the per-alignment table

```bash
# into your active AMR++ environment, or any python3 env you point the pipeline at
pip install pysam duckdb
```

**R** (for `plot_sweep_dropoff.R`):

- `data.table`
- `ggplot2`
- `scales`
- `patchwork` for the two-panel heat maps

```bash
# from an R session, or Rscript -e '...'
install.packages(c("data.table", "ggplot2", "scales", "patchwork"))
```

> If any of these are missing, the step that needs them will fail. The Python step is required to produce any output. The R step only affects the figures, so without R you still get the combined CSV matrices. If `patchwork` alone is missing, the plots are still written, but the count and percentage heat maps come out as separate files rather than paired panels.

Make sure `Rscript` is on your `PATH`. The pipeline calls it through the `$RSCRIPT` environment variable, so set `RSCRIPT` in your config or profile if your Rscript lives somewhere non-standard.

---

## What the sweep does

For each BAM, the pipeline:

1. Streams every primary alignment once and records, per alignment, how much of the read aligned, how well it matched, and which reference positions it covered.
2. Re-counts the resistome at every combination in the threshold grid, using the same fragment-level, group-aware, edge-aware counting as the main pipeline.
3. Writes per-sample result CSVs.

Then across all samples:

4. Combines the per-sample CSVs into single matrices with `combine_sweep_results.py`.
5. Produces drop-off plots and summary tables with `plot_sweep_dropoff.R`.

Because the BAM is parsed only once and the whole grid is evaluated from the resulting intermediate table, the full grid costs roughly the same as a single counting run.

A sample that ends up with no usable alignments, for example one where every alignment was to a gene requiring SNP confirmation and `sweep_exclude_snp` removed it, is handled gracefully. The sweep writes header-only outputs for that sample, prints a warning naming the filters responsible, and the run continues.

---

## Running the sweep

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline bam_coverage_sweep \
    --bam_files "/path/to/Alignment/BAM_files/Standard/*_alignment_sorted.bam" \
    --output Filtering_sensitivity_analysis
```

- `--bam_files` is a glob matching the BAMs to sweep. Quote it so the shell does not expand it.
- `--output` is where results are written.
- On a cluster, swap `-profile local` for your SLURM profile. See [Running with SLURM](Running_with_SLURM.md).

Example with the SLURM profile:

```bash
nextflow run main_AMR++.nf \
    --pipeline bam_coverage_sweep \
    --bam_files "AMR++_output/Alignment/BAM_files/Standard/*" \
    --output Filtering_sensitivity_analysis \
    -profile local_slurm
```

---

## The three filters

AMR++ applies three independent alignment filters. The sweep explores the same three, so a threshold you pick from a sweep can be applied directly to a counting run.

**All values are proportions from 0 to 1.** A value of `0` turns a filter off.

| Filter | What it asks |
|---|---|
| gene fraction | What proportion of the **gene** was covered by retained alignments? |
| query coverage | What proportion of the **read** aligned? |
| identity | What proportion of the **aligned portion** matched the reference? |

These are independent rather than substitutes for one another. A short fragment that aligns perfectly scores high on identity and low on query coverage. A full-length read carrying many mismatches scores the reverse.

### Two ways to measure query coverage

Query coverage has two definitions, and you choose between them rather than applying both. The threshold value means the same thing either way. Only the metric it is compared against changes.

| Parameter | Calculation | Behavior |
|---|---|---|
| `sweep_query_coverage` | `aligned_length / read_length` | Counts mismatches inside the aligned region as covered |
| `sweep_match_qcov` | `(aligned_length - NM) / read_length` | Counts only genuine matches |

Because the edit distance is never negative, the match-based value can never exceed the standard one, so a given threshold is always at least as strict under `sweep_match_qcov`.

In the counting step the same choice is made with the `match_qcov` parameter, which is `"Y"` or `"N"` rather than a list of values.

### The NM tag requirement

`sweep_match_qcov` and `sweep_identity` are both computed from the NM tag, which records each alignment's edit distance. NM is an optional SAM field, so it is not always present.

With the filter off, alignments lacking NM still pass. At any threshold above 0 they are **excluded**, because their identity cannot be verified. BWA-MEM emits NM by default, so BAMs produced by the standard AMR++ alignment step are fine, but BAMs rewritten by other tools may not carry it.

If you need to add it:

```bash
samtools calmd -b input.bam reference.fasta > output.nm.bam
samtools index output.nm.bam
```

---

## Choosing what to sweep

The sweep always tests a **two-dimensional grid**: gene fraction against one read-level filter.

**Axis 1 is always gene fraction.** `sweep_gene_fraction` is swept in every run and appears in every figure.

**Axis 2 is one of the three read-level filters.** Set your choice to a list of values and leave the others at `0`.

The relevant `params.config` block:

```groovy
    /* AXIS 1: always swept */
    sweep_gene_fraction   = "0,0.1,0.25,0.5,0.8"

    /* AXIS 2: sweep ONE of these three; leave the others at "0",
     * or give one a single value to apply it as a fixed filter */
    sweep_query_coverage  = "0,0.5,0.6,0.7,0.8,0.9,0.95"
    sweep_match_qcov      = "0"
    sweep_identity        = "0"

    /* Applied at every point in the grid regardless of what is swept */
    sweep_edge_aware_qcov = "Y"
    sweep_exclude_snp     = "N"
```

`sweep_query_coverage` and `sweep_match_qcov` measure the same property, so sweeping both is meaningless. Pick whichever definition of query coverage you want to filter on. `sweep_identity` measures something different and is a genuine alternative axis.

### A single value is a fixed filter

A single **non-zero** value is a filter applied at every point in the grid rather than a swept axis. This is how you hold one filter constant while exploring another:

```groovy
// sweep gene fraction against query coverage, with identity pinned at 0.9
sweep_gene_fraction   = "0,0.1,0.25,0.5,0.8"
sweep_query_coverage  = "0,0.5,0.6,0.7,0.8,0.9,0.95"
sweep_identity        = "0.9"
```

The figures plot gene fraction against query coverage, and the plot subtitle records that identity was fixed at 0.9 so the constraint stays visible.

### Overriding on the command line

```bash
nextflow run main_AMR++.nf --pipeline bam_coverage_sweep \
    --bam_files "AMR++_output/Alignment/BAM_files/Standard/*" \
    --sweep_gene_fraction '0,0.5,0.9' \
    --sweep_query_coverage '0,0.6,0.7,0.8,0.85,0.9,0.95' \
    --output Filtering_sensitivity_analysis
```

Quote the comma-separated values. Most shells handle them unquoted, but quoting removes any ambiguity, particularly inside SLURM submission scripts.

The plotting script reads the gene-fraction grid from the results themselves, so changing the grid needs no corresponding edit anywhere else.

### Grid size

Every combination is evaluated, so the grid grows multiplicatively. Rows written to `combined_results.csv` per sample:

| Configuration | Grid | Rows per sample |
|---|---|---|
| Gene fraction against query coverage (default) | 5 x 7 | 35 |
| Gene fraction against identity | 5 x 5 | 25 |
| A finer grid, 8 gene fraction by 10 query coverage | 8 x 10 | 80 |

`combined_gene_detail.csv` multiplies these by the number of detected genes, so it grows quickly. Runtime scales far better than file size, since the BAM is parsed once, but keep an eye on disk usage for large sample sets.

---

## Outputs

Under `--output`, in `CoverageSweep/`:

- `PerSample/` contains one set of CSVs per BAM: `<prefix>_results.csv`, `_gene_detail.csv`, `_redundancy.csv` and `_length_quantiles.csv`.
- `Combined/` contains the merged matrices across all samples: `combined_results.csv`, `combined_gene_detail.csv` and `combined_length_quantiles.csv`, each with a `sample_id` column.
- `Summary/` contains the drop-off plots as PNG files and the summary tables as CSV.

---

## Interpreting the plots

**Drop-off curves** (`dropoff_<level>_count_vs_primary.png` and `_pct_vs_primary.png`) show, for each annotation level, how the count falls as the read-level threshold tightens. The bold line is the mean across samples and each thin grey line is one sample. One line is drawn per gene-fraction value, so both axes appear together.

**Heat maps** (`dropoff_<level>_heatmap.png` and `reads_retained_combined_heatmap.png`) show the full two-dimensional grid as a pair of panels: the absolute feature count on the left and the percentage of baseline retained on the right. The percentage panel uses a fixed 0 to 100 colour scale, so colour carries the same meaning across every level and every figure.

**Reads retained and lost** figures show what fraction of fragment hits survive each threshold relative to the no-filter baseline.

Expect the annotation levels to behave differently. Gene accession counts usually fall fastest, because accession-level detection is the most sensitive to database redundancy, while group, mechanism and class counts hold up longer. A threshold that removes many accessions but few groups is generally preferable to the reverse, since aggregation to the group level or above is what we recommend for statistical comparison in any case.

The sweep shows what a threshold costs on your data. It cannot tell you which threshold is correct, because that depends on the question you are asking.

---

## Applying the selected filter thresholds

Once you have chosen thresholds, re-run counting on the same BAM files with `bam_resistome`. The parameter names match the sweep parameters without the `sweep_` prefix, and use the same 0 to 1 proportion scale.

| Sweep parameter | Counting parameter |
|---|---|
| `sweep_gene_fraction` | `--min_gene_fraction` |
| `sweep_query_coverage` | `--min_query_coverage` |
| `sweep_identity` | `--min_identity` |
| `sweep_match_qcov` | `--match_qcov Y` together with `--min_query_coverage` |

A run using a gene fraction of 0.5 and a query coverage of 0.8:

```bash
nextflow run main_AMR++.nf \
    --pipeline bam_resistome \
    --bam_files "AMR++_output/Alignment/BAM_files/Standard/*bam" \
    --min_gene_fraction 0.5 \
    --min_query_coverage 0.8 \
    --output AMR++_gf50_qcov80 \
    --prefix AMR_gf50_qcov80 \
    -profile local_slurm
```

If your sweep used `sweep_match_qcov` rather than `sweep_query_coverage`, set `match_qcov` when applying so the counting step computes query coverage the same way the sweep did:

```bash
nextflow run main_AMR++.nf \
    --pipeline bam_resistome \
    --bam_files "AMR++_output/Alignment/BAM_files/Standard/*bam" \
    --min_gene_fraction 0.5 \
    --min_query_coverage 0.8 \
    --match_qcov Y \
    --output AMR++_gf50_matchqcov80 \
    --prefix AMR_gf50_matchqcov80 \
    -profile local_slurm
```

This produces a new count matrix under `Results/` reflecting the selected filters, ready for aggregation and statistical analysis.

> ### Keep a record of which filters produced which matrix
>
> **The count matrix filename does not encode the filter settings.** Every run writes
> `<prefix>_analytic_matrix.csv`, which by default is `AMR_analytic_matrix.csv` no matter which
> thresholds were applied. Two runs with different filters produce two files with identical
> names, and if they share an output directory the second silently overwrites the first.
> Nothing inside the file records how it was generated.
>
> Three habits prevent this from becoming a problem.
>
> **Give each run its own `--output` directory.** The simplest safeguard, and it also keeps the
> per-sample intermediate files separated.
>
> **Encode the thresholds in `--prefix`.** The prefix becomes part of the matrix filename, so
> `--prefix AMR_gf50_qcov80` yields `AMR_gf50_qcov80_analytic_matrix.csv`, which stays
> self-describing even after the file is moved or shared. Both examples above do this.
>
> **Save the command.** Keep the exact invocation alongside the results, or move the settings
> into a params file you can archive:
>
> ```bash
> nextflow run main_AMR++.nf -params-file my_thresholds.yaml --pipeline bam_resistome ...
> ```
