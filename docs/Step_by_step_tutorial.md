# Paired-End Analysis: Step-by-Step Tutorial

This tutorial walks through a complete paired-end AMR++ analysis, running each step independently. This approach is recommended for large datasets because it lets you delete the `work/` directory between steps to manage storage, and makes it easier to troubleshoot individual steps.

> **Before you start:** Make sure you have run the demo at least once from a login node (`nextflow run main_AMR++.nf -profile local --pipeline demo`). This installs the SNP confirmation software and verifies your environment. See [Getting Started](GettingStarted.md) for details.

## Table of Contents

- [Setup](#setup)
  - [Load the Conda Environment](#load-the-conda-environment)
  - [Confirm Tools Are Available](#confirm-tools-are-available)
- [Step 1: Quality Assessment (eval_qc)](#step-1-quality-assessment-eval_qc)
- [Step 2: Quality Trimming (trim_qc)](#step-2-quality-trimming-trim_qc)
- [Step 2.5 (Optional): Read Deduplication (dedup)](#step-25-optional-read-deduplication-dedup)
- [Step 3: Host Read Removal (rm_host)](#step-3-host-read-removal-rm_host)
- [Step 4: Resistome Analysis (resistome)](#step-4-resistome-analysis-resistome)
- [Step 5 (Optional): Microbiome Analysis (kraken)](#step-5-optional-microbiome-analysis-kraken)
- [Managing the Work Directory](#managing-the-work-directory)
- [Effect of Key Parameters](#effect-of-key-parameters)

---

## Setup

### Load the Conda Environment

AMR++ requires a set of bioinformatic tools that are bundled into a conda environment. You only need to create the environment once.

```bash
# If you haven't installed the conda environment yet
conda env create -f envs/AMR++_env.yaml

# Activate it before each session
conda activate AMR++_env
```

If your HPC uses modules to provide conda, load it first:

```bash
module load Anaconda3/2024.02-1
conda activate AMR++_env
```

### Confirm Tools Are Available

```bash
nextflow -version
bwa 2>&1 | head -3
samtools --version | head -1
trimmomatic -version
```

If any of these fail, review the [installation document](installation.md) or check that your conda environment is activated.

---

## Step 1: Quality Assessment (eval_qc)

**What it does:** Runs FastQC on all input reads and generates a MultiQC summary report so you can assess read quality before trimming.

**Input:** Raw paired-end FASTQ files  
**Output:** `AMR++_results/QC_analysis/FastQC/` and `AMR++_results/QC_analysis/MultiQC_stats/multiqc_report.html`

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline eval_qc \
    --reads "data/raw/*_R{1,2}.fastq.gz" \
    --output AMR++_results
```

**Explore the results:**
```bash
ls AMR++_results/QC_analysis/FastQC/
# Download and open multiqc_report.html in a browser to review quality across all samples
```

**When done:** Review the MultiQC report. Look for adapter contamination, low-quality tails, or unusual GC content — these inform your trimming parameters in Step 2.

```bash
# Safe to delete the work directory after confirming output
rm -rf work/
```

---

## Step 2: Quality Trimming (trim_qc)

**What it does:** Runs Trimmomatic to remove adapter sequences and low-quality bases from reads.

**Input:** Raw paired-end FASTQ files  
**Output:**
- `AMR++_results/QC_trimming/Paired/` — trimmed paired reads (`*.1P.fastq.gz`, `*.2P.fastq.gz`)
- `AMR++_results/QC_trimming/Unpaired/` — reads whose pair was discarded
- `AMR++_results/Results/Stats/trimmomatic.stats` — read counts before/after trimming

**Default trimming parameters** (adjust in `params.config` or on the command line as needed):

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
    --output AMR++_results
```

**Inspect results:**
```bash
ls AMR++_results/QC_trimming/Paired/
cat AMR++_results/Results/Stats/trimmomatic.stats
```

```bash
rm -rf work/
```

---

## Step 2.5 (Optional): Read Deduplication (dedup)

**What it does:** Uses seqkit to remove exact duplicate reads by sequence. We recommend this step for **target-enriched sequencing data** (e.g., probe-based capture), where PCR duplicates are more prevalent. For standard shotgun metagenomics, this step is typically skipped.

**Input:** QC-trimmed paired reads  
**Output:** `AMR++_results/Deduped_reads/` — deduplicated reads (`*_R1.dedup.fastq.gz`, `*_R2.dedup.fastq.gz`)

> **Note:** Deduplication runs on R1 and R2 independently. A read pair is removed if R1 is an exact duplicate of another R1. See the [Analysis Recommendations](Analysis_recommendations.md) for more context.

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline dedup \
    --reads "AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz" \
    --output AMR++_results
```

**Inspect results:**
```bash
ls AMR++_results/Deduped_reads/
```

```bash
rm -rf work/
```

---

## Step 3: Host Read Removal (rm_host)

**What it does:** Aligns reads against the host genome using BWA. Reads that map to the host are discarded; unmapped reads are kept for downstream analysis.

**Input:** QC-trimmed (or deduped) paired reads  
**Output:**
- `AMR++_results/HostRemoval/NonHostFastq/` — non-host reads (`*.non.host.R1.fastq.gz`, `*.non.host.R2.fastq.gz`)
- `AMR++_results/Results/Stats/host.removal.stats` — read counts before/after removal

**Parameters to set:**

| Parameter | Description |
|-----------|-------------|
| `--host` | Path to your host genome FASTA. For bovine on Grace HPRC: `/scratch/group/big_scratch/SHARED_resources/host_genome/GCF_002263795.3_ARS-UCD2.0_genomic.fna` |
| `--reads` | Point to trimmed reads (or deduped reads if you ran Step 2.5) |

```bash
# Using QC-trimmed reads
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline rm_host \
    --reads "AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --output AMR++_results

# OR using deduped reads (if you ran Step 2.5)
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline rm_host \
    --reads "AMR++_results/Deduped_reads/*R{1,2}.dedup.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --output AMR++_results
```

**Inspect results:**
```bash
ls AMR++_results/HostRemoval/NonHostFastq/
cat AMR++_results/Results/Stats/host.removal.stats
```

```bash
rm -rf work/
```

---

## Step 4: Resistome Analysis (resistome)

**What it does:** Aligns reads to the MEGARes antimicrobial resistance gene database using BWA. Generates count matrices at multiple annotation levels (Type, Class, Mechanism, Group, Gene), with optional SNP confirmation and deduplication of alignments.

**Input:** Non-host FASTQ files from Step 3  
**Output:**
- `AMR++_results/Results/AMR_analytic_matrix.csv` — standard count matrix
- `AMR++_results/Results/SNPconfirmed_AMR_analytic_matrix.csv` — SNP-confirmed counts (if `--snp Y`)
- `AMR++_results/Results/dedup_AMR_analytic_matrix.csv` — deduplicated counts (if `--deduped Y`)
- `AMR++_results/ResistomeAnalysis/ResistomeCounts/` — per-sample count files

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threshold` | 0 | Gene fraction threshold. We recommend 0 and aggregating to Group level for analysis |
| `--snp` | N | Set to `Y` to enable SNP confirmation for relevant genes |
| `--deduped` | N | Set to `Y` to also output deduplicated alignment counts |

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline resistome \
    --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" \
    --snp Y \
    --deduped Y \
    --output AMR++_results
```

**Inspect results:**
```bash
ls AMR++_results/Results/
head AMR++_results/Results/AMR_analytic_matrix.csv
wc -l AMR++_results/Results/AM*  # compare row counts across output matrices
```

```bash
rm -rf work/
```

---

## Step 5 (Optional): Microbiome Analysis (kraken)

**What it does:** Classifies non-host reads taxonomically using Kraken2 and generates a wide-format count matrix across all samples.

**Input:** Non-host FASTQ files from Step 3  
**Output:**
- `AMR++_results/MicrobiomeAnalysis/Kraken/` — per-sample raw and report files
- `AMR++_results/Results/kraken_analytic_matrix.csv` — aggregated taxonomic matrix

**Parameters to set:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kraken_db` | null | Path to your Kraken2 database directory |
| `--kraken_confidence` | 0.0 | Classification confidence threshold (0–1). Higher values = fewer but more confident classifications |

> **Memory note:** Kraken2 loads the entire database into RAM by default. For a large database (e.g., core_nt, ~200+ GB), ensure you are running with a SLURM `xlarge` label or equivalent memory allocation. See [Running with SLURM](Running_with_SLURM.md).

```bash
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline kraken \
    --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" \
    --kraken_db /path/to/kraken_db \
    --output AMR++_results
```

**Inspect results:**
```bash
head AMR++_results/Results/kraken_analytic_matrix.conf_0.0.csv
head AMR++_results/Results/unclassifieds_kraken_analytic_matrix.conf_0.0.csv
```

**Effect of confidence threshold:**
```bash
# Compare classification rates at different confidence levels
nextflow run main_AMR++.nf \
    -profile local \
    --pipeline kraken \
    --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" \
    --kraken_db /path/to/kraken_db \
    --kraken_confidence 0.5 \
    --output AMR++_results

head AMR++_results/Results/unclassifieds_kraken_analytic_matrix.conf_0.5.csv
```

```bash
rm -rf work/
```

---

## Managing the Work Directory

As mentioned above, the `work/` directory stores all intermediate files and enables `-resume`. Here is a practical strategy for managing it:

```bash
# Check how large the work directory has grown
du -sh work/

# Delete it safely after a step completes successfully
rm -rf work/

# Redirect to a scratch drive with more space
nextflow run main_AMR++.nf -profile local --pipeline trim_qc \
    --reads "data/raw/*_R{1,2}.fastq.gz" \
    -w /scratch/$USER/nf_work_trim
```

If a step fails partway through, use `-resume` to pick up where it left off (the work directory must still exist):

```bash
nextflow run main_AMR++.nf -profile local --pipeline rm_host \
    --reads "AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --output AMR++_results \
    -resume
```

---

## Effect of Key Parameters

### Gene Fraction Threshold

The `--threshold` parameter controls the minimum fraction of a gene that must be covered by alignments for it to be counted. The default is 0 (count all alignments). We recommend keeping this at 0 and instead aggregating your count matrix to the **Group** level during statistical analysis to account for cross-mapping between closely related genes.

```bash
# Compare results at different thresholds
nextflow run main_AMR++.nf --pipeline resistome \
    --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" \
    --threshold 80 --output AMR++_thresh80_results

wc -l AMR++_results/Results/AMR_analytic_matrix.csv
wc -l AMR++_thresh80_results/Results/AMR_analytic_matrix.csv
```

### Kraken Confidence

```bash
# A confidence of 1.0 will classify almost nothing — useful to see extreme
head AMR++_results/Results/unclassifieds_kraken_analytic_matrix.conf_1.0.csv
```

See [Analysis Recommendations](Analysis_recommendations.md) for guidance on interpreting these results.
