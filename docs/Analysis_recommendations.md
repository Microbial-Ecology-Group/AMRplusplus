# Analysis Recommendations

This document contains tips and recommendations for running AMR++ and interpreting your results. These suggestions are based on our experience analyzing resistome and microbiome data across various sequencing approaches and sample types.

## Table of Contents

- [AMR++ Pipeline Run](#amr-pipeline-run)
  - [Profile Use](#profile-use)
  - [Deduplicated Counts](#deduplicated-counts)
- [Resistome Statistics](#resistome-statistics)
  - [Gene Fraction Threshold](#gene-fraction-threshold)
  - [Alignment Counting](#alignment-counting)
  - [Rarefaction Analysis](#rarefaction-analysis)
- [Microbiome Statistics](#microbiome-statistics)



# AMR++ Pipeline Run

## Profile Use

If you have SLURM as a job scheduler, we recommend using a profile ending in `_slurm` as this can facilitate running a lot of samples, especially if you are using a large Kraken database. More details on how to [run AMR++ with SLURM](Running_with_SLURM.md).

- In our case, we use NCBI's NT database which can require as much as 250 GB of memory for classification.
- If you use a `local` profile but need to submit with `sbatch`, the memory allocation will be split across the number of ongoing processes (controlled by `maxForks`). This can cause a process to fail if multiple high-memory processes are running at once.
- By using a SLURM profile, AMR++ can automatically submit individual jobs with smaller resource requirements when possible and submit short-duration, high-memory processes when needed for Kraken2. (Note: we could use Kraken2's `--memory-mapping` flag to reduce memory usage, but we find this takes too long for large samples.)

## Deduplicated Counts

We typically recommend deduplicated alignment counts for target-enriched sequencing samples.

---

# Resistome Statistics

## Gene Fraction Threshold

As of AMR++ v4, we have relaxed the gene fraction `--threshold` from 80% to 0%.

The gene fraction threshold was originally applied with the goal of reducing false positive calls for individual gene accessions. Instead, we now recommend including all gene accessions in the final count matrix and then aggregating counts to the AMR **Group** level. Genes with high sequence homology (e.g., beta-lactam resistance genes) are most likely to have misclassifications at the individual gene level. Aggregating to the Group level accounts for this because Group annotations were created based on reference sequence homology.

> **Example:** We might not be able to differentiate between the individual genes "acc-1" and "acc-2", but they are both in the "ACC" group.

Any characterization of individual genes should be performed with alternative bioinformatic methods.

## Alignment Counting

As of AMR++ v4, we only count 1 alignment hit per read. Previously we counted all alignments (including secondary and supplemental), but we now recommend only keeping primary alignments. You can still change this by modifying the `--samtools_flag` parameter.

## Rarefaction Analysis

The rarefaction analysis is helpful for determining whether you achieved adequate sequencing depth for the resistome. This is different from the `rarefaction()` function in R, because with AMR++ we also incorporate the probability of a read *not* being classified as an ARG gene. In comparison, using `rarefaction()` in R with your count matrix will only subset from reads already classified as ARGs.

---

# Microbiome Statistics

## Microbiome counts
For the paired-end and merged-read analyses, we opted to only count each pair of reads as a single hit:

- **Paired-end reads** — We use Kraken2's `--paired` flag.
- **Merged reads** — We count classifications for merged reads normally, but include the `--interleaved` flag for the unmerged reads.
- **Single-end reads** — We count each classification normally.

## Kraken2 database
Based on our experience, we find that using the most comprehensive database possible provides the best results (maximizing sensitivity/reducing false-positive calls). We currently use the `core_nt` database as [maintained by Ben Langmead](https://benlangmead.github.io/aws-indexes/k2). 