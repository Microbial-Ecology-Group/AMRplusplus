# Choosing the right pipeline

AMR++ analyzes data by combining workflows that take sequencing reads through various bioinformatic software. We recommend our standard AMR++ pipeline as a comprehensive way to start from raw sequencing reads, perform QC assessment, host DNA removal, and resistome analysis with MEGARes. However, users might only want to run portions of the pipeline or have more control over their computing needs. 

Using the `--pipeline` parameter, users can select from complete pipelines or individual subworkflows.

## Table of Contents
- [Quick Reference](#quick-reference)
- [Pipeline Selection Flowchart](#pipeline-selection-flowchart)
- [Main Pipelines (Paired-End Reads)](#main-pipelines-paired-end-reads)
  - [Paired-End Subworkflows](#paired-end-subworkflows)
- [Merged Read Pipelines](#merged-read-pipelines)
  - [Merged Read Subworkflows](#merged-read-subworkflows)
- [Single-End Pipeline](#single-end-pipeline)
  - [Single-End Subworkflows](#single-end-subworkflows)
- [BAM File Workflows](#bam-file-workflows)
- [Complete Example Command](#complete-example-command)

## Quick Reference

| Data Type | Complete Pipeline | Description |
|-----------|------------------|-------------|
| Paired-end | `standard_AMR` | Full analysis with host removal |
| Paired-end | `fast_AMR` | Fast analysis, no host removal |
| Paired-end | `standard_AMR_wKraken` | Full analysis + microbiome |
| Merged reads | `merged_AMR` | For FLASH-merged paired reads |
| Merged reads | `merged_AMR_wKraken` | Merged reads + microbiome |
| Single-end | `se_AMR` | Single-end reads analysis |
| Single-end | `se_AMR_wKraken` | Single-end reads + microbiome |
| BAM files | `bam_resistome` | From pre-aligned BAMs |

## Pipeline Selection Flowchart

```
What type of sequencing data do you have?
│
├── Paired-end reads (R1 + R2)
│   ├── Need host removal? 
│   │   ├── Yes → standard_AMR or standard_AMR_wKraken
│   │   └── No  → fast_AMR
│   └── Already merged with FLASH?
│       └── Yes → merged_AMR or merged_AMR_wKraken
│
├── Single-end reads
│   └── se_AMR
│
└── Pre-aligned BAM files
    └── bam_resistome or bam_resistome_counts
```


---

## Main Pipelines (Paired-End Reads)

These pipelines process standard paired-end Illumina reads using the `--reads` parameter.

### `--pipeline demo` (default)
Simple demonstration using included test data. Runs `fast_AMR` workflow.
```bash
nextflow run main_AMR++.nf -profile local --pipeline demo
```

### `--pipeline standard_AMR`
**Recommended for most analyses.** Complete pipeline with host removal.
- **Steps:** QC trimming → Host DNA removal → Resistome alignment → Resistome results
```bash
nextflow run main_AMR++.nf -profile local --pipeline standard_AMR \
    --reads "path/to/reads/*_R{1,2}.fastq.gz" \
    --host /path/to/host_genome.fasta
```

### `--pipeline fast_AMR`
Streamlined pipeline that skips host removal. Use when host contamination is minimal or not a concern.
- **Steps:** QC trimming → Resistome alignment → Resistome results
```bash
nextflow run main_AMR++.nf -profile local --pipeline fast_AMR \
    --reads "path/to/reads/*_R{1,2}.fastq.gz"
```

### `--pipeline standard_AMR_wKraken`
Full pipeline with added microbiome analysis using Kraken2.
- **Steps:** 
  - QC trimming → Host DNA removal → Resistome alignment → Resistome results
  - Non-host reads → Kraken2 taxonomic classification
- **Note:** Requires a Kraken database. The minikraken_8GB_202003 database (~8GB) will be downloaded automatically, or specify your own with `--kraken_db`.
```bash
nextflow run main_AMR++.nf -profile local --pipeline standard_AMR_wKraken \
    --reads "path/to/reads/*_R{1,2}.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --kraken_db /path/to/kraken_db/
```

## Paired-End Subworkflows

Individual analysis steps that can be run independently on paired-end data.

| Subworkflow | Description | Key Parameters |
|-------------|-------------|----------------|
| `eval_qc` | FastQC quality assessment | `--reads` |
| `trim_qc` | Trimmomatic adapter/quality trimming | `--reads` |
| `rm_host` | BWA host read removal | `--reads`, `--host` |
| `resistome` | Full resistome analysis | `--reads`, `--amr`, `--annotation` |
| `align` | BWA alignment to MEGARes only | `--reads`, `--amr` |
| `kraken` | Kraken2 taxonomic classification | `--reads`, `--kraken_db` |
| `qiime2` | QIIME 2 16S rRNA analysis | `--reads`, `--dada2_db` |

### Examples

```bash
# Quality assessment only
nextflow run main_AMR++.nf -profile local --pipeline eval_qc \
    --reads "path/to/reads/*_R{1,2}.fastq.gz"

# Trimming only
nextflow run main_AMR++.nf -profile local --pipeline trim_qc \
    --reads "path/to/reads/*_R{1,2}.fastq.gz"

# Host removal only
nextflow run main_AMR++.nf -profile local --pipeline rm_host \
    --reads "path/to/reads/*_R{1,2}.fastq.gz" \
    --host /path/to/host_genome.fasta

# Resistome analysis (alignment + counting + rarefaction)
nextflow run main_AMR++.nf -profile local --pipeline resistome \
    --reads "path/to/reads/*_R{1,2}.fastq.gz"

# Alignment only (generates BAM files)
nextflow run main_AMR++.nf -profile local --pipeline align \
    --reads "path/to/reads/*_R{1,2}.fastq.gz"

# Kraken classification
nextflow run main_AMR++.nf -profile local --pipeline kraken \
    --reads "path/to/reads/*_R{1,2}.fastq.gz" \
    --kraken_db /path/to/kraken_db/

# QIIME 2 16S analysis
nextflow run main_AMR++.nf -profile local --pipeline qiime2 \
    --reads "path/to/reads/*_R{1,2}.fastq.gz" \
    --dada2_db /path/to/dada2_db/
```


---

## Merged Read Pipelines

For paired-end reads that have been merged using FLASH or similar tools. Use the `--merged_reads` parameter to specify input files.

### `--pipeline merged_AMR`
Standard pipeline for FLASH-merged paired-end reads. Processes both merged and unmerged read fractions.
- **Steps:** QC trimming → Host DNA removal → Resistome alignment → Resistome results
```bash
nextflow run main_AMR++.nf -profile local --pipeline merged_AMR \
    --merged_reads "path/to/merged/*_{merged,unmerged}.fastq.gz" \
    --host /path/to/host_genome.fasta
```

### `--pipeline merged_AMR_wKraken`
Merged read pipeline with Kraken2 microbiome analysis.
```bash
nextflow run main_AMR++.nf -profile local --pipeline merged_AMR_wKraken \
    --merged_reads "path/to/merged/*_{merged,unmerged}.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --kraken_db /path/to/kraken_db/
```

## Merged Read Subworkflows

Individual analysis steps for merged paired-end reads.

| Subworkflow | Description | Key Parameters |
|-------------|-------------|----------------|
| `merge_reads` | Merge paired-end reads with FLASH | `--reads` |
| `merged_rm_host` | Host removal for merged reads | `--merged_reads`, `--host` |
| `merged_resistome` | Resistome analysis for merged reads | `--merged_reads`, `--amr`, `--annotation` |
| `merged_kraken` | Kraken2 for merged reads | `--merged_reads`, `--kraken_db` |

### Examples

```bash
# Merge paired-end reads with FLASH
nextflow run main_AMR++.nf -profile local --pipeline merge_reads \
    --reads "path/to/reads/*_R{1,2}.fastq.gz"

# Host removal on merged reads
nextflow run main_AMR++.nf -profile local --pipeline merged_rm_host \
    --merged_reads "path/to/merged/*.extendedFrags.fastq.gz" \
    --host /path/to/host_genome.fasta

# Resistome analysis on merged reads
nextflow run main_AMR++.nf -profile local --pipeline merged_resistome \
    --merged_reads "path/to/merged/*_{merged,unmerged}.fastq.gz"
```

---

## Single-End Pipeline

For single-end sequencing data. Uses the `--reads` parameter with single files.

### `--pipeline se_AMR`
Complete AMR++ pipeline for single-end reads with Kraken analysis.
- **Steps:** QC trimming → Host DNA removal → Resistome alignment → Resistome results
```bash
nextflow run main_AMR++.nf -profile local --pipeline se_AMR \
    --reads "path/to/reads/*.fastq.gz" \
    --host /path/to/host_genome.fasta
```

### `--pipeline se_AMR_wKraken`
Complete AMR++ pipeline for single-end reads with Kraken analysis.
- **Steps:** QC trimming → Host DNA removal → Resistome alignment → Resistome results → Kraken2
```bash
nextflow run main_AMR++.nf -profile local --pipeline se_AMR \
    --reads "path/to/reads/*.fastq.gz" \
    --host /path/to/host_genome.fasta
```

## Single-End Subworkflows

Individual analysis steps for single-end data.

| Subworkflow | Description | Key Parameters |
|-------------|-------------|----------------|
| `se_eval_qc` | FastQC for single-end reads | `--reads` |
| `se_trim_qc` | Trimmomatic for single-end reads | `--reads` |
| `se_rm_host` | Host removal for single-end reads | `--reads`, `--host` |
| `se_resistome` | Resistome analysis for single-end reads | `--reads`, `--amr`, `--annotation` |
| `se_kraken` | Kraken2 for single-end reads | `--reads`, `--kraken_db` |

### Examples

```bash
# Single-end QC
nextflow run main_AMR++.nf -profile local --pipeline se_eval_qc \
    --reads "path/to/reads/*.fastq.gz"

# Single-end resistome
nextflow run main_AMR++.nf -profile local --pipeline se_resistome \
    --reads "path/to/reads/*.fastq.gz"
```


---

## BAM File Workflows

For analyzing pre-aligned BAM files (e.g., from a previous AMR++ alignment run).

| Workflow | Description | Key Parameters |
|----------|-------------|----------------|
| `bam_resistome` | Full resistome analysis from BAM files | `--bam_files`, `--amr`, `--annotation` |
| `bam_resistome_counts` | Generate count matrices from BAM files | `--bam_files`, `--amr`, `--annotation` |

### Examples

```bash
# Resistome analysis from existing BAM files
nextflow run main_AMR++.nf -profile local --pipeline bam_resistome \
    --bam_files "path/to/alignments/*.bam"

# Count matrix generation from BAM files
nextflow run main_AMR++.nf -profile local --pipeline bam_resistome_counts \
    --bam_files "path/to/alignments/*.bam"
```

---

## Complete Example Command

In the following example, we run the standard AMR++ workflow with SNP verification and deduplicated counts:

```bash
nextflow run main_AMR++.nf -profile local --pipeline standard_AMR \
    --reads "path/to/your/reads/*_R{1,2}.fastq.gz" \
    --host /path/to/host_genome.fasta \
    --snp Y \
    --deduped Y
```

This will generate:
- `AMR_analytic_matrix.csv` - Standard count matrix
- `SNPconfirmed_AMR_analytic_matrix.csv` - SNP-confirmed counts
- `Deduped_AMR_analytic_matrix.csv` - Deduplicated counts
- `Deduped_SNPconfirmed_AMR_analytic_matrix.csv` - Deduplicated + SNP-confirmed counts

---

