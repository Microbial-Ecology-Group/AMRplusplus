# Coverage Threshold Sweep: Step-by-Step Tutorial

This tutorial covers the `bam_coverage_sweep` pipeline, which takes a set of pre-aligned BAM files and evaluates the resistome across a grid of **gene-fraction** and **query-coverage** thresholds. Instead of committing to a single cutoff, the sweep shows you how the number of detected genes (and classes, mechanisms, and groups) drops off as you tighten each threshold — so you can choose sensible cutoffs and see how sensitive each sample is to them. It produces combined result matrices and a set of drop-off plots.

> **Before you start:** This pipeline runs on BAM files you have already produced with AMR++ (for example, the `*_alignment_sorted.bam` files under `Alignment/BAM_files/Standard/`). It also requires a few extra tools that are **not** bundled in the default AMR++ conda environment — see [Additional tools required](#additional-tools-required) below.

## Table of Contents

- [Additional tools required](#additional-tools-required)
- [What the sweep does](#what-the-sweep-does)
- [Running the sweep](#running-the-sweep)
- [Choosing the threshold values](#choosing-the-threshold-values)
- [Outputs](#outputs)
- [Interpreting the plots](#interpreting-the-plots)

---

## Additional tools required

The sweep uses a Python analysis step and an R plotting step whose dependencies are **not
included in the AMR++ conda environment**, because adding them (especially the R plotting
stack) would substantially increase an already-large environment. Install them separately
before running the sweep.

**Python (for `coverage_threshold_sweep.py` and `combine_sweep_results.py`):**
- `pysam` — reading BAM alignments
- `duckdb` — fast on-disk aggregation of the per-alignment table

```bash
# into your active AMR++ environment (or any python3 env you point the pipeline at)
pip install pysam duckdb
```

**R (for `plot_sweep_dropoff.R`):**
- `data.table`
- `ggplot2`
- `scales`

```bash
# from an R session, or Rscript -e '...'
install.packages(c("data.table", "ggplot2", "scales"))
```

> If any of these are missing, the sweep step that needs them will fail. The Python step is
> required to produce any output; the R step only affects the plots — if R isn't available you
> will still get the combined CSV matrices, just not the figures.

Make sure `Rscript` is on your `PATH` (the pipeline calls it via the `$RSCRIPT` environment
variable; set `RSCRIPT` in your config/profile if your Rscript is in a non-standard location).

---

## What the sweep does

For each BAM, the pipeline:
1. Streams every primary alignment and records, per gene, its gene-fraction (breadth of the
   reference covered) and query-coverage (fraction of the read aligned).
2. Re-counts the resistome at every combination of the gene-fraction and query-coverage
   thresholds in the sweep grid — using the same fragment-level, group-aware, edge-aware
   counting as the main pipeline.
3. Writes per-sample result CSVs.

Then across all samples it:
4. Combines the per-sample CSVs into single matrices (`combine_sweep_results.py`).
5. Produces drop-off plots and summary tables (`plot_sweep_dropoff.R`).

---

## Running the sweep

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline bam_coverage_sweep \
    --bam_files "/path/to/Alignment/BAM_files/Standard/*_alignment_sorted.bam" \
    --output MSS_sweep
```

- `--bam_files` is a glob matching the BAMs to sweep. Quote it so the shell doesn't expand it.
- `--output` is where results are written.
- On a cluster, swap `-profile local` for your SLURM profile (e.g. `-profile local_slurm`) —
  see [Running with SLURM](Running_with_SLURM.md).

Example with the SLURM profile:
```bash
nextflow run main_AMR++.nf \
    --pipeline bam_coverage_sweep \
    --bam_files "/scratch/$USER/MSS_output/Alignment/BAM_files/Standard/*" \
    --output MSS_sweep \
    -profile local_slurm
```

---

## Choosing the threshold values

The grids are defined as defaults in `bin/coverage_threshold_sweep.py` and can be overridden
on the command line. Both gene-fraction and query-coverage are given as **proportions (0–1)**:

- `--gene-fraction-sweep` default: `0, 0.1, 0.25, 0.5, 0.8`
- `--query-coverage-sweep` default: `0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95`
- `--identity-sweep` default: `0` (off); note that identity, if used, is on a **0–100 percent**
  scale, not a proportion.

The default grid evaluates every combination (7 × 5 = 35 threshold pairs) in a single pass.

> If you change the gene-fraction grid, update `GF_SWEEP` near the top of
> `bin/plot_sweep_dropoff.R` to match, so the plots use the same values.

---

## Outputs

Under `--output`, in `CoverageSweep/`:

- `PerSample/` — one set of CSVs per BAM (`<prefix>_results.csv`, `_gene_detail.csv`,
  `_length_quantiles.csv`).
- `Combined/` — the merged matrices across all samples (`combined_results.csv`,
  `combined_gene_detail.csv`, `combined_length_quantiles.csv`), each with a `sample_id` column.
- `Summary/` (figures directory) — the drop-off plots (PNG) and summary tables (CSV).

---

## Interpreting the plots

- **Drop-off curves** (`dropoff_<level>_count_vs_primary`, `_pct_vs_primary`) show, for each
  annotation level (gene, class, mechanism, group), how the count falls as the primary
  threshold tightens. The bold line is the mean across samples; each thin grey line is one
  sample.
- **Heatmaps** (`dropoff_<level>_heatmap`, `reads_retained_combined_heatmap_*`) show the joint
  effect of both thresholds at once.
- **Reads retained / lost** figures show what fraction of fragment-hits survive each threshold
  relative to the no-filter baseline.

Use these to pick thresholds where the resistome is stable (a plateau) rather than in a steep
drop-off region, and to spot samples that are unusually sensitive to filtering. Once you've
chosen cutoffs, apply them in the main pipeline via `--min_gene_fraction` and
`--min_query_coverage` (same 0–1 proportion scale).