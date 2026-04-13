# Single-End Analysis: Step-by-Step Tutorial

This tutorial walks through a complete AMR++ analysis for single-end sequencing data, running each step independently. The single-end pipeline mirrors the paired-end workflow but uses SE-specific subworkflows prefixed with `se_`.

> **Before you start:** Make sure you have run the demo at least once from a login node (`nextflow run main_AMR++.nf -profile local --pipeline demo`). This installs the SNP confirmation software and verifies your environment. See [Getting Started](GettingStarted.md) for details.

## Table of Contents

- [Setup](#setup)
  - [Load the Conda Environment](#load-the-conda-environment)
  - [Single-End vs Paired-End Input](#single-end-vs-paired-end-input)
- [Step 1: Quality Assessment (eval_qc)](#step-1-quality-assessment-eval_qc)
- [Step 2: Quality Trimming (se_trim_qc)](#step-2-quality-trimming-se_trim_qc)
- [Step 2.5 (Optional): Read Deduplication (se_dedup)](#step-25-optional-read-deduplication-se_dedup)
- [Step 3: Host Read Removal (se_rm_host)](#step-3-host-read-removal-se_rm_host)
- [Step 4: Resistome Analysis (se_resistome)](#step-4-resistome-analysis-se_resistome)
- [Step 5 (Optional): Microbiome Analysis (se_kraken)](#step-5-optional-microbiome-analysis-se_kraken)
- [Managing the Work Directory](#managing-the-work-directory)

---

## Setup

### Load the Conda Environment

```bash
# If you haven't installed the conda environment yet
conda env create -f envs/AMR++_env.yaml

# Activate before each session
conda activate AMR++_env
```

On HPC systems using modules:
```bash
module load Anaconda3/2024.02-1
conda activate AMR++_env
```

Confirm tools are available:
```bash
nextflow -version
bwa 2>&1 | head -3
samtools --version | head -1
```

### Single-End vs Paired-End Input

The key difference in running the single-end pipeline is how you specify `--reads`. For single-end data, each sample is a single FASTQ file, so the glob pattern does not use the `{1,2}` syntax:

```bash
# Paired-end (matches R1 and R2 pairs)
--reads "data/raw/*_R{1,2}.fastq.gz"

# Single-end (matches individual files)
--reads "data/raw/*.fastq.gz"
```

AMR++ automatically detects the read type from the glob pattern and routes to the correct pipeline.

---

## Step 1: Quality Assessment (eval_qc)

**What it does:** Runs FastQC on all input reads and generates a MultiQC summary report.

**Input:** Raw single-end FASTQ files  
**Output:** `SE_AMR++_results/QC_analysis/FastQC/` and `SE_AMR++_results/QC_analysis/MultiQC_stats/multiqc_report.html`

> Note: `eval_qc` works the same for both single-end and paired-end data — you do not need a separate `se_eval_qc` pipeline for this step.

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline eval_qc \
    --reads "data/raw/*.fastq.gz" \
    --output SE_AMR++_results
```

**Explore results:**
```bash
ls SE_AMR++_results/QC_analysis/FastQC/
# Download and open multiqc_report.html in a browser
```

```bash
rm -rf work/
```

---

## Step 2: Quality Trimming (se_trim_qc)

**What it does:** Runs Trimmomatic in SE mode to remove adapters and low-quality bases.

**Input:** Raw single-end FASTQ files  
**Output:**
- `SE_AMR++_results/QC_trimming/Single/` — trimmed reads (`*.trimmed.fastq.gz`)
- `SE_AMR++_results/Results/Stats/trimmomatic.stats` — per-sample read counts

**Default trimming parameters** (modify in `params.config` or on the command line):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--adapters` | nextera.fa | Adapter sequences to remove |
| `--leading` | 3 | Minimum quality for leading bases |
| `--trailing` | 3 | Minimum quality for trailing bases |
| `--slidingwindow` | 4:15 | Window size : required average quality |
| `--minlen` | 36 | Minimum read length after trimming |

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_trim_qc \
    --reads "data/raw/*.fastq.gz" \
    --output SE_AMR++_results
```

**Inspect results:**
```bash
ls SE_AMR++_results/QC_trimming/Single/
cat SE_AMR++_results/Results/Stats/trimmomatic.stats
```

```bash
rm -rf work/
```

---

## Step 2.5 (Optional): Read Deduplication (se_dedup)

**What it does:** Uses seqkit to remove exact duplicate reads by sequence. Recommended for **target-enriched sequencing data** (probe-based capture) where PCR duplicates are common.

**Input:** QC-trimmed single-end reads  
**Output:** `SE_AMR++_results/Deduped_reads/` — deduplicated reads (`*.dedup.fastq.gz`)

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_dedup \
    --reads "SE_AMR++_results/QC_trimming/Single/*.trimmed.fastq.gz" \
    --output SE_AMR++_results
```

**Inspect results:**
```bash
ls SE_AMR++_results/Deduped_reads/
```

```bash
rm -rf work/
```

---

## Step 3: Host Read Removal (se_rm_host)

**What it does:** Aligns single-end reads against the host genome using BWA. Unmapped reads are extracted and saved as non-host FASTQ files for downstream analysis.

**Input:** Trimmed (or deduped) single-end reads  
**Output:**
- `SE_AMR++_results/HostRemoval/NonHostFastq/` — non-host reads (`*.non.host.fastq.gz`)
- `SE_AMR++_results/Results/Stats/host.removal.stats` — read counts

**Parameters to set:**

| Parameter | Description |
|-----------|-------------|
| `--host` | Path to your host genome FASTA. For bovine on Grace HPRC: `/scratch/group/big_scratch/SHARED_resources/host_genome/GCF_002263795.3_ARS-UCD2.0_genomic.fna` |
| `--reads` | Point to trimmed reads (or deduped reads if you ran Step 2.5) |

```bash
# Using QC-trimmed reads
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_rm_host \
    --reads "SE_AMR++_results/QC_trimming/Single/*.trimmed.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --output SE_AMR++_results

# OR using deduped reads (if you ran Step 2.5)
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_rm_host \
    --reads "SE_AMR++_results/Deduped_reads/*.dedup.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --output SE_AMR++_results
```

**Inspect results:**
```bash
ls SE_AMR++_results/HostRemoval/NonHostFastq/
cat SE_AMR++_results/Results/Stats/host.removal.stats
```

```bash
rm -rf work/
```

---

## Step 4: Resistome Analysis (se_resistome)

**What it does:** Aligns non-host reads to the MEGARes database and generates count matrices at all annotation levels, with optional SNP confirmation and deduplicated counts.

**Input:** Non-host FASTQ files from Step 3  
**Output:**
- `SE_AMR++_results/Results/AMR_analytic_matrix.csv`
- `SE_AMR++_results/Results/SNPconfirmed_AMR_analytic_matrix.csv` (if `--snp Y`)
- `SE_AMR++_results/Results/dedup_AMR_analytic_matrix.csv` (if `--deduped Y`)
- `SE_AMR++_results/ResistomeAnalysis/ResistomeCounts/` — per-sample count files

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threshold` | 0 | Gene fraction threshold. Keep at 0 and aggregate to Group level for analysis |
| `--snp` | N | Set to `Y` to enable SNP confirmation |
| `--deduped` | N | Set to `Y` to output deduplicated alignment counts |

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_resistome \
    --reads "SE_AMR++_results/HostRemoval/NonHostFastq/*.non.host.fastq.gz" \
    --snp Y \
    --deduped Y \
    --output SE_AMR++_results
```

**Inspect results:**
```bash
ls SE_AMR++_results/Results/
head SE_AMR++_results/Results/AMR_analytic_matrix.csv
wc -l SE_AMR++_results/Results/AM*
```

```bash
rm -rf work/
```

---

## Step 5 (Optional): Microbiome Analysis (se_kraken)

**What it does:** Classifies non-host reads taxonomically using Kraken2 and generates a wide-format count matrix.

**Input:** Non-host FASTQ files from Step 3  
**Output:**
- `SE_AMR++_results/MicrobiomeAnalysis/Kraken/` — per-sample reports
- `SE_AMR++_results/Results/kraken_analytic_matrix.conf_${kraken_confidence}.csv`

> **Memory note:** Kraken2 loads the entire database into RAM. For large databases, make sure your job has sufficient memory. See [Running with SLURM](Running_with_SLURM.md) for memory allocation details.

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_kraken \
    --reads "SE_AMR++_results/HostRemoval/NonHostFastq/*.non.host.fastq.gz" \
    --kraken_db /path/to/kraken_db \
    --output SE_AMR++_results
```

**Inspect results:**
```bash
head SE_AMR++_results/Results/kraken_analytic_matrix.conf_0.0.csv
head SE_AMR++_results/Results/unclassifieds_kraken_analytic_matrix.conf_0.0.csv
```

```bash
rm -rf work/
```

---

## Managing the Work Directory

The `work/` directory enables `-resume` but grows large quickly. Delete it after each step completes successfully.

```bash
# Check size
du -sh work/

# Delete after confirming output looks correct
rm -rf work/

# Redirect to a scratch drive with more space
nextflow run main_AMR++.nf -profile local --pipeline se_rm_host \
    --reads "SE_AMR++_results/QC_trimming/Single/*.trimmed.fastq.gz" \
    --host /path/to/host_genome.fasta \
    -w /scratch/$USER/nf_work_host_removal
```

To resume a failed run (work directory must still exist):
```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline se_resistome \
    --reads "SE_AMR++_results/HostRemoval/NonHostFastq/*.non.host.fastq.gz" \
    --output SE_AMR++_results \
    -resume
```
