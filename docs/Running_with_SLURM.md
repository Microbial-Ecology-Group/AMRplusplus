# AMR++ SLURM Resource Labels

## Overview

When running AMR++ with a SLURM-enabled profile (any profile ending in `_slurm`, such as `local_slurm`, `conda_slurm`, or `singularity_slurm`), each process in the pipeline is assigned a resource label that determines how much memory, time, and CPUs it requests from the SLURM scheduler. This allows individual jobs to request only the resources they need rather than reserving a single large allocation for the entire pipeline.

To use this setup, submit the Nextflow command itself inside an `sbatch` script that requests **minimal resources and a long walltime** — for example, 1 CPU, 4 GB of memory, and 24–48 hours. This lightweight "driver" job stays running on the cluster while Nextflow handles submitting and monitoring each pipeline step as a separate SLURM job. The `queueSize` parameter (set to 10 in `local_slurm.config`) controls how many jobs Nextflow is allowed to have queued or running at the same time, preventing the pipeline from flooding the scheduler.

Example sbatch wrapper:
```bash
#!/bin/bash
#SBATCH --job-name=AMRplusplus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00

nextflow run main_AMR++.nf -profile conda_slurm --pipeline standard_AMR \
    --reads 'data/raw/*_R{1,2}.fastq.gz'
```

Most resource tiers use 20 CPUs by default and include automatic retry logic — if a job fails with an out-of-memory exit code (137–140), it retries once with 1.5× the memory and time.

---

## Resource Tiers

The resource labels are defined in `config/local_slurm.config`. Each label specifies a combination of memory, walltime, and retry behavior.

### `nano` — 512 MB / 5 min

Minimal resources for tasks that download files or copy small dependencies.

| Process | Module |
|---------|--------|
| `build_dependencies` | Resistome/resistome.nf |

---

### `micro` — 10 GB / 30 min

Lightweight tasks: indexing, QC reporting, plotting, and result aggregation.

| Process | Module |
|---------|--------|
| `index` | Alignment/bwa.nf |
| `HostRemovalStats` | Alignment/bwa.nf |
| `fastqc` | Fastqc/fastqc.nf |
| `multiqc` | Fastqc/fastqc.nf |
| `dlkraken` | Microbiome/kraken2.nf |
| `krakenresults` | Microbiome/kraken2.nf |
| `runbracken` | Microbiome/kraken2.nf |
| `plotrarefaction` | Resistome/resistome.nf |
| `snpresults` | Resistome/resistome.nf |

---

### `micro_long` — 10 GB / 1.5 hr

Same memory as `micro` but with extended walltime for I/O-heavy trimming of large FASTQ files.

| Process | Module |
|---------|--------|
| `runqc` (paired-end Trimmomatic) | Trimming/trimmomatic.nf |

---

### `small` — 20 GB / 2 hr

Moderate tasks: alignment, read merging, deduplication, and summary statistics.

| Process | Module |
|---------|--------|
| `bwa_align` | Alignment/bwa.nf |
| `bwa_align_se` | Alignment/bwa.nf |
| `bwa_merged_align` | Alignment/bwa.nf |
| `samtools_dedup_se` | Alignment/bwa.nf |
| `samtools_merge_bams` | Alignment/bwa.nf |
| `MergeReadsFlash` | QC/merge.nf |
| `SeqkitReadCounts` | QC/merge.nf |
| `Qiime2Import` | Microbiome/qiime2.nf |
| `resistomeresults` | Resistome/resistome.nf |
| `runrarefaction` | Resistome/resistome.nf |
| `runqc_se` (SE Trimmomatic) | Trimming/trimmomatic.nf |
| `QCstats` | Trimming/trimmomatic.nf |
| `QCstats_SE` | Trimming/trimmomatic.nf |

---

### `medium` — 24 GB / 6 hr

Resource-intensive tasks: host removal, resistome counting, and QIIME 2 analysis steps.

| Process | Module |
|---------|--------|
| `bwa_rm_contaminant_fq` | Alignment/bwa.nf |
| `bwa_rm_contaminant_merged_fq` | Alignment/bwa.nf |
| `bwa_rm_contaminant_se` | Alignment/bwa.nf |
| `runresistome` | Resistome/resistome.nf |
| `Qiime2Dada2` | Microbiome/qiime2.nf |
| `Qiime2Classify` | Microbiome/qiime2.nf |
| `Qiime2Filter` | Microbiome/qiime2.nf |
| `Qiime2Tree` | Microbiome/qiime2.nf |
| `Qiime2Export` | Microbiome/qiime2.nf |

---

### `large` — 32 GB / 24 hr

Currently not assigned to any process, but available for future use or manual override.

---

### `xlarge` — 256 GB / 1.5 hr

High-memory, short-duration tasks. This is the default tier for Kraken2 when the database is loaded entirely into RAM (no `--memory-mapping` flag).

| Process | Module |
|---------|--------|
| `runkraken`* | Microbiome/kraken2.nf |
| `runkraken_merged`* | Microbiome/kraken2.nf |
| `runkraken_se`* | Microbiome/kraken2.nf |

*\*These processes use a dynamic label — they switch to `large` (32 GB / 24 hr) if `--memory-mapping` is included in `params.kraken_options`, since memory-mapped mode uses much less RAM but takes longer.*

---

### `snp_ignore` — 20 GB / 2 hr (errors ignored)

Special tier for SNP verification, which may fail for samples without SNP-associated genes. Failures are silently ignored so they don't halt the pipeline.

| Process | Module |
|---------|--------|
| `runsnp` | Resistome/resistome.nf |

---

## Additional Labels

Two labels referenced in some config files but not currently assigned to any process:

- **`large_short`** — 50 GB / 2 hr
- **`large`** — 32 GB / 24 hr

These are available for future processes or can be used with Nextflow's `-process.label` CLI override.

---

## Customizing Resources

To adjust resources for your cluster, edit `config/local_slurm.config`. For example, to reduce the CPUs for all tiers:

```groovy
withLabel: small {
    cpus = 10    // changed from 20
    memory = { 20.GB * (task.attempt > 1 ? 1.5 : 1) }
    time = { 2.h * (task.attempt > 1 ? 1.5 : 1) }
}
```

You can also override specific labels from the command line without editing the config:

```bash
nextflow run main_AMR++.nf -profile conda_slurm \
    --pipeline standard_AMR \
    -process.memory '64 GB' \
    -process.time '12h'
```
