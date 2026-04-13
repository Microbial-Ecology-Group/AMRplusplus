# Merged-Read Analysis: Step-by-Step Tutorial

This tutorial covers the AMR++ workflow for paired-end reads that are merged using FLASH before analysis. The merged-read approach is useful for amplicon or short-insert libraries where a significant proportion of R1 and R2 reads overlap and can be joined into a single longer fragment. Both the merged (overlapping) and unmerged (non-overlapping) reads are carried through the pipeline separately and combined at the results stage.

> **Before you start:** Make sure you have run the demo at least once from a login node (`nextflow run main_AMR++.nf -profile local --pipeline demo`). This installs the SNP confirmation software and verifies your environment. Your conda environment must also include FLASH and seqKit — if you installed the AMR++ environment from `envs/AMR++_env.yaml`, these are already included. See [Getting Started](GettingStarted.md) for details.

## Table of Contents

- [Setup](#setup)
  - [Load the Conda Environment](#load-the-conda-environment)
  - [How Merged-Read Analysis Works](#how-merged-read-analysis-works)
- [Step 1: Quality Assessment (eval_qc)](#step-1-quality-assessment-eval_qc)
- [Step 2: Quality Trimming (trim_qc)](#step-2-quality-trimming-trim_qc)
- [Step 2.5 (Optional): Read Deduplication (dedup)](#step-25-optional-read-deduplication-dedup)
- [Step 3: Merge Reads with FLASH (merge_reads)](#step-3-merge-reads-with-flash-merge_reads)
- [Step 4: Host Read Removal (merged_rm_host)](#step-4-host-read-removal-merged_rm_host)
- [Step 5: Resistome Analysis (merged_resistome)](#step-5-resistome-analysis-merged_resistome)
- [Step 6 (Optional): Microbiome Analysis (merged_kraken)](#step-6-optional-microbiome-analysis-merged_kraken)
- [Managing the Work Directory](#managing-the-work-directory)

---

## Setup

### Load the Conda Environment

```bash
# Create the environment if you haven't already
conda env create -f envs/AMR++_env.yaml

# Activate before each session
conda activate AMR++_env
```

On HPC systems using modules:
```bash
module load Anaconda3/2024.02-1
conda activate AMR++_env
```

Confirm that FLASH and seqKit are available:
```bash
flash --version
seqkit version
nextflow -version
```

### How Merged-Read Analysis Works

After trimming, paired-end reads are passed to FLASH, which attempts to overlap R1 and R2 reads. This produces two output files per sample:

| File suffix | Description |
|-------------|-------------|
| `.extendedFrags.fastq.gz` | Reads where R1 and R2 overlapped and were merged into a single longer fragment |
| `.notCombined.fastq.gz` | Reads where R1 and R2 did not overlap and remain as separate SE-like reads |

From Step 4 onward, both the merged and unmerged read fractions are processed in parallel and their results are combined. All downstream steps use `--merged_reads` instead of `--reads` to point to these two-file-per-sample outputs.

> **Important:** When specifying `--merged_reads`, always use **single quotes** (`'`), not backticks or double quotes, around glob patterns that contain curly braces:
> ```bash
> --merged_reads 'path/to/files/*.{extendedFrags,notCombined}.fastq.gz'
> ```

---

## Step 1: Quality Assessment (eval_qc)

**What it does:** Runs FastQC and generates a MultiQC report on the raw paired-end reads before any processing.

**Input:** Raw paired-end FASTQ files  
**Output:** `Merged_AMR++_results/QC_analysis/`

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline eval_qc \
    --reads "data/raw/*_R{1,2}.fastq.gz" \
    --output Merged_AMR++_results
```

```bash
# Review the MultiQC report, then clean up
ls Merged_AMR++_results/QC_analysis/
rm -rf work/
```

---

## Step 2: Quality Trimming (trim_qc)

**What it does:** Runs Trimmomatic in PE mode to remove adapters and low-quality bases before merging.

**Input:** Raw paired-end FASTQ files  
**Output:**
- `Merged_AMR++_results/QC_trimming/Paired/` — trimmed paired reads (`*.1P.fastq.gz`, `*.2P.fastq.gz`)
- `Merged_AMR++_results/Results/Stats/trimmomatic.stats`

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
    --pipeline trim_qc \
    --reads "data/raw/*_R{1,2}.fastq.gz" \
    --output Merged_AMR++_results
```

```bash
ls Merged_AMR++_results/QC_trimming/Paired/
rm -rf work/
```

---

## Step 2.5 (Optional): Read Deduplication (dedup)

**What it does:** Uses seqkit to remove exact duplicate reads by sequence. Recommended for **target-enriched sequencing data**. This step is run on the trimmed paired-end reads *before* merging.

**Input:** QC-trimmed paired reads  
**Output:** `Merged_AMR++_results/Deduped_reads/`

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline dedup \
    --reads "Merged_AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz" \
    --output Merged_AMR++_results
```

```bash
ls Merged_AMR++_results/Deduped_reads/
rm -rf work/
```

---

## Step 3: Merge Reads with FLASH (merge_reads)

**What it does:** Attempts to overlap and merge R1 and R2 reads using FLASH. Produces two output files per sample: merged (extended fragments) and unmerged (not combined).

**Input:** QC-trimmed paired reads (or deduped reads if you ran Step 2.5)  
**Output:** `Merged_AMR++_results/Flash_reads/` — one `*.extendedFrags.fastq.gz` and one `*.notCombined.fastq.gz` per sample

```bash
# Using QC-trimmed reads
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline merge_reads \
    --reads "Merged_AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz" \
    --output Merged_AMR++_results

# OR using deduped reads (if you ran Step 2.5)
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline merge_reads \
    --reads "Merged_AMR++_results/Deduped_reads/*R{1,2}.dedup.fastq.gz" \
    --output Merged_AMR++_results
```

**Inspect results:**
```bash
ls Merged_AMR++_results/Flash_reads/
# You should see pairs of files per sample:
#   sample.extendedFrags.fastq.gz
#   sample.notCombined.fastq.gz
```

```bash
rm -rf work/
```

---

## Step 4: Host Read Removal (merged_rm_host)

**What it does:** Aligns both the merged and unmerged read fractions against the host genome. Non-host reads are extracted from each fraction separately.

**Input:** FLASH output files from Step 3 (use `--merged_reads`, not `--reads`)  
**Output:** `Merged_AMR++_results/HostRemoval/NonHostFastq/` — one `*.merged.non.host.fastq.gz` and one `*.unmerged.non.host.fastq.gz` per sample

**Parameters to set:**

| Parameter | Description |
|-----------|-------------|
| `--host` | Path to your host genome FASTA. For bovine on Grace HPRC: `/scratch/group/big_scratch/SHARED_resources/host_genome/GCF_002263795.3_ARS-UCD2.0_genomic.fna` |
| `--merged_reads` | Glob pointing to FLASH output files — **use single quotes** |

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline merged_rm_host \
    --merged_reads 'Merged_AMR++_results/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz' \
    --host /path/to/host_genome.fasta \
    --output Merged_AMR++_results
```

**Inspect results:**
```bash
ls Merged_AMR++_results/HostRemoval/NonHostFastq/
# Expected output per sample:
#   sample.merged.non.host.fastq.gz
#   sample.unmerged.non.host.fastq.gz
```

```bash
rm -rf work/
```

---

## Step 5: Resistome Analysis (merged_resistome)

**What it does:** Aligns both the merged and unmerged non-host fractions to MEGARes and generates combined count matrices with optional SNP confirmation.

**Input:** Non-host FASTQ files from Step 4 (use `--merged_reads`, not `--reads`)  
**Output:**
- `Merged_AMR++_results/Results/AMR_analytic_matrix.csv`
- `Merged_AMR++_results/Results/SNPconfirmed_AMR_analytic_matrix.csv` (if `--snp Y`)
- `Merged_AMR++_results/Results/dedup_AMR_analytic_matrix.csv` (if `--deduped Y`)

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threshold` | 0 | Gene fraction threshold. Keep at 0 and aggregate to Group level |
| `--snp` | N | Set to `Y` to enable SNP confirmation |
| `--deduped` | N | Set to `Y` to output deduplicated alignment counts (different than read deduplication) |

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline merged_resistome \
    --merged_reads 'Merged_AMR++_results/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz' \
    --snp Y \
    --deduped Y \
    --output Merged_AMR++_results
```

**Inspect results:**
```bash
ls Merged_AMR++_results/Results/
head Merged_AMR++_results/Results/AMR_analytic_matrix.csv
wc -l Merged_AMR++_results/Results/AM*
```

```bash
rm -rf work/
```

---

## Step 6 (Optional): Microbiome Analysis (merged_kraken)

**What it does:** Classifies both merged and unmerged non-host reads using Kraken2. Merged reads are classified normally; unmerged reads use the `--interleaved` flag so each pair counts as a single classification. Results are combined into a single matrix.

**Input:** Non-host FASTQ files from Step 4 (use `--merged_reads`)  
**Output:**
- `Merged_AMR++_results/MicrobiomeAnalysis/Kraken/` — per-sample reports for merged and unmerged
- `Merged_AMR++_results/Results/kraken_analytic_matrix.conf_${kraken_confidence}.csv`

> **Memory note:** Kraken2 requires enough RAM to load the full database. For large databases (e.g., core_nt), use a SLURM `xlarge` allocation. See [Running with SLURM](Running_with_SLURM.md).

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline merged_kraken \
    --merged_reads 'Merged_AMR++_results/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz' \
    --kraken_db /path/to/kraken_db \
    --output Merged_AMR++_results
```

**Inspect results:**
```bash
head Merged_AMR++_results/Results/kraken_analytic_matrix.conf_0.0.csv
head Merged_AMR++_results/Results/unclassifieds_kraken_analytic_matrix.conf_0.0.csv
```

```bash
rm -rf work/
```

---

## Managing the Work Directory

The `work/` directory enables `-resume` but can grow to many gigabytes between steps. Delete it after each step completes successfully.

```bash
# Check size
du -sh work/

# Delete after confirming output is correct
rm -rf work/

# Redirect to a scratch drive
nextflow run main_AMR++.nf -profile local --pipeline merge_reads \
    --reads "Merged_AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz" \
    -w /scratch/$USER/nf_work_merge
```

To resume a failed run (work directory must still exist):
```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline merged_rm_host \
    --merged_reads 'Merged_AMR++_results/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz' \
    --host /path/to/host_genome.fasta \
    --output Merged_AMR++_results \
    -resume
```
