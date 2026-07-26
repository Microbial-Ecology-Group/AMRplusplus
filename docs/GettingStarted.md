# Getting Started with AMR++

This document will help you get AMR++ running on your system and understand the key concepts you'll encounter throughout your analysis. Before diving into a full tutorial, it's worth reading this page to understand how the pipeline is structured and how to get the most out of it.

## Table of Contents

- [Quick Start](#quick-start)
- [Key Concepts](#key-concepts)
  - [The Work Directory](#the-work-directory)
  - [Running the Pipeline in Steps](#running-the-pipeline-in-steps)
  - [Running the Demo First](#running-the-demo-first)
  - [Conda vs Local Profiles](#conda-vs-local-profiles)
  - [Running with SLURM](#running-with-slurm)
- [Tutorials](#tutorials)
- [Other Helpful Documents](#other-helpful-documents)

---

## Quick Start

If you just want to get something running, here are the minimum steps:

```bash
# 1. Clone the AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git
cd AMRplusplus

# 2. Create and activate the conda environment
conda env create -f envs/AMR++_env.yaml
conda activate AMR++_env

# 3. Run the demo to confirm everything works and install SNP software
nextflow run main_AMR++.nf -profile local --pipeline demo
```

If the demo completes without errors, you're ready to analyze your data. Continue reading below to understand how to configure the pipeline for your needs.

---

## Key Concepts

### The Work Directory

Every time you run AMR++, Nextflow creates a `work/` directory in your current folder. This directory stores all intermediate files for every task — it's what allows Nextflow to `-resume` a failed run from where it left off without recomputing everything from scratch.

**Important things to know about the work directory:**

- It can grow very large, often tens to hundreds of gigabytes for a full dataset
- It is safe to delete it once a pipeline step completes successfully, *as long as you don't need to `-resume` that run*
- You can redirect it to a different location (e.g., a scratch drive) using the `-w` flag

```bash
# Redirect the work directory to a scratch drive with more space
nextflow run main_AMR++.nf --pipeline trim_qc -w /scratch/$USER/nf_work
```

**Recommended practice when running in steps:** delete the work directory between steps to avoid filling your storage. For example:

```bash
nextflow run main_AMR++.nf --pipeline trim_qc --output my_results
# Confirm output looks correct, then:
rm -rf work/
```

---

### Running the Pipeline in Steps

AMR++ supports running individual analysis steps independently using the `--pipeline` flag. This is recommended for large datasets because:

1. **Storage management** — you can delete the `work/` directory between steps
2. **Troubleshooting** — it's easier to identify and fix problems one step at a time
3. **Flexibility** — you can swap in your own outputs at any stage (e.g., use pre-trimmed reads)
4. **Compute efficiency** — on HPC systems, you only request resources for the specific step you need

The typical order of steps for a paired-end analysis is:

```
eval_qc → trim_qc → (dedup) → rm_host → resistome → (kraken)
```

For single-end data:
```
eval_qc → se_trim_qc → (se_dedup) → se_rm_host → se_resistome → (se_kraken)
```

For merged-read data:
```
eval_qc → trim_qc → (dedup) → merge_reads → merged_rm_host → merged_resistome → (merged_kraken)
```

See [Choosing the Right Pipeline](choosing_pipeline.md) for a full list of pipeline options and when to use each.

---

### Running the Demo First

**Always run the demo before your first real analysis.** The demo serves two purposes:

1. Confirms that your environment and all software dependencies are working correctly
2. Downloads and installs the [SNP confirmation software](https://github.com/Isabella136/AmrPlusPlus_SNP) into `bin/AmrPlusPlus_SNP/`

The SNP software is cloned from GitHub during the `build_dependencies` step. On HPC systems, compute nodes often do not have internet access, which means if you skip the demo and go straight to a full pipeline run on a compute node, the SNP software download will fail.

```bash
# Run the demo once from a login node (which has internet access)
nextflow run main_AMR++.nf -profile local --pipeline demo
```

After the demo completes, the SNP software will be cached locally in `bin/` and subsequent pipeline runs on compute nodes will find it there without needing internet access.

If you submit jobs via `sbatch` and skip the demo, you can alternatively add a WebProxy module to your sbatch script:

```bash
# Add to your sbatch script before the nextflow command
module load WebProxy
```

---

### Conda vs Local Profiles

AMR++ uses the `-profile` flag to determine how it finds the required bioinformatic tools. The two most common options are `conda` and `local`, and understanding the difference matters for how you set up your jobs.

**`-profile conda`**

Nextflow manages the conda environment itself — it will create and activate it automatically. This is convenient but means Nextflow needs conda available at runtime.

```bash
nextflow run main_AMR++.nf -profile conda --pipeline demo
```

**`-profile local`**

AMR++ assumes all required tools are already available in your current `$PATH`. This is the profile to use if you have already activated the conda environment yourself before running the pipeline.

```bash
# Activate the environment yourself first
conda activate AMR++_env

# Then run with -profile local — tools are already in your PATH
nextflow run main_AMR++.nf -profile local --pipeline demo
```

**When does this matter?** On HPC systems with SLURM, you typically activate the conda environment in your `sbatch` script before calling nextflow, then use `-profile local_slurm`. This gives you more control over which environment is active and avoids Nextflow trying to manage conda in a job context. Example sbatch script:

```bash
#!/bin/bash
#SBATCH --job-name=AMRplusplus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00

# Activate conda environment before calling nextflow
module load Anaconda3/2024.02-1
conda activate AMR++_env

# Use local_slurm — tools are in PATH, SLURM handles job submission
nextflow run main_AMR++.nf -profile local_slurm \
    --pipeline standard_AMR \
    --reads 'data/raw/*_R{1,2}.fastq.gz'
```

---

### Running with SLURM

If your institution uses SLURM as a job scheduler, using a `_slurm` profile is strongly recommended, especially for large datasets or when using a large Kraken database.

**How it works:**

When you use a `_slurm` profile (e.g., `local_slurm`, `conda_slurm`, `singularity_slurm`), you submit a lightweight "driver" job to SLURM — typically requesting just 1 CPU, 4 GB of memory, and a long walltime (24–48 hours). The Nextflow driver process stays running and automatically submits each individual pipeline step as its own separate SLURM job with the appropriate resource requests.

This means:
- High-memory steps like Kraken2 (which may need 256 GB) are submitted with that allocation automatically
- Low-memory steps like indexing or result aggregation use minimal resources
- You don't need to request the maximum memory for all steps upfront

```bash
#!/bin/bash
#SBATCH --job-name=AMRplusplus_driver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00

module load Anaconda3/2024.02-1
conda activate AMR++_env

nextflow run main_AMR++.nf -profile local_slurm \
    --pipeline standard_AMR \
    --reads 'data/raw/*_R{1,2}.fastq.gz' \
    --host /path/to/host_genome.fasta
```

**Controlling parallelism:** The `queueSize` parameter in the SLURM config (default: 10) controls how many jobs Nextflow submits to the scheduler simultaneously. You can raise this if your cluster allows more concurrent jobs.

**Memory retry:** If a job fails with an out-of-memory error (exit codes 137–140), AMR++ automatically retries it with 1.5× the memory allocation — up to 3 attempts.

For full details on resource tiers and how to customize them, see [Running with SLURM](Running_with_SLURM.md).

---

**Running without SLURM (local mode):**

If you are on a single workstation or a login node with enough resources, you can run the entire pipeline locally:

```bash
nextflow run main_AMR++.nf -profile local \
    --pipeline standard_AMR \
    --reads 'data/raw/*_R{1,2}.fastq.gz'
```

In this mode, all processes run sequentially or with limited parallelism on the same machine. The `maxForks` setting in the config controls how many processes run simultaneously.

---

## Tutorials

Step-by-step instructions for each read type:

- **[Paired-end reads](Step_by_step_tutorial.md)** — standard Illumina R1 + R2 reads
- **[Single-end reads](SingleEnd_read_tutorial.md)** — single-file FASTQ input
- **[Merged reads](Merged_read_tutorial.md)** — reads merged with FLASH before analysis

---

## Other Helpful Documents

| Document | Description |
|----------|-------------|
| [Installation](installation.md) | How to install AMR++ and its dependencies |
| [Choosing a Pipeline](choosing_pipeline.md) | Full list of `--pipeline` options and when to use each |
| [Configuration](configuration.md) | All parameters, how to modify `params.config`, profiles |
| [Running with SLURM](Running_with_SLURM.md) | Resource labels, sbatch setup, memory tuning |
| [Output](output.md) | Directory structure and description of all output files |
| [Analysis Recommendations](Analysis_recommendations.md) | Tips on gene fraction thresholds, counting, rarefaction |
| [Dependencies](dependencies.md) | Full list of software tools used |
| [FAQs & Troubleshooting](FAQs.md) | Common errors and how to fix them |
| [CHANGELOG](CHANGELOG.md) | Version history and updates |
