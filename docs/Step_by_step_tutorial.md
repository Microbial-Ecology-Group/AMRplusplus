
# Table of contents
* [Evaluate your server](#evaluate-your-server)
* [Run AMR++ demo](#run-first-demo)
* [Evaluate sample QC](#run-the-eval_qc-pipeline)
* [Run QC trimming](#run-trim_qc-pipeline)
* [Run host removal](#run-host-removal-pipeline)
* [Run resistome](#run-resistome-pipeline)
* [Run microbiome](#run-kraken-pipeline)

# Evaluate your server

Let's explore your server and identify what tools are available.

```
ls

pwd

cd workshop/data

ls

# Check which tools are installed
bwa --help
```

We don't have all the tools installed, but we can use "Anaconda" to load an environment with all necessary tools. 


```
module load Anaconda3/2024.02-1

source activate /home/training/conda_envs/AMR++_env
```


# Run first Demo

Let's download the code for [AMR++ from github](https://github.com/Microbial-Ecology-Group/AMRplusplus/tree/master) using the "git" command.


```
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

cd AMRplusplus/

ls
```


Next, we can run the demo with this simple command:

```
nextflow run main_AMR++.nf
```

## How did AMR++ know what parameters to run:

AMR++ has default parameters that are listed in the `params.config` file. We can then either change the values in the file directly, or add those flags to the command you're using directly. We'll practice that below. 

For now, let's look at the defaults:
``` 
less params.config 
```


## Running AMR++ in pieces

This just ran a demonstration of AMR++ using the included test data. This is extremely useful for small datasets or if you have access to large computing clusters. However, it is often computationally prohibitive to run the entire pipeline so we'll now learn how to run each component individually. 


# Run the eval_qc pipeline

Input files for `--pipeline eval_qc`:
* sample reads (`--reads`) (by default `--reads` = ${baseDir}/data/raw/*_R{1,2}.fastq.gz)

Output files:
* FastQC results for each sample
* MultiQC report with aggregated results

## Run first pipeline, eval_qc

For this subworkflow, we'll just specify the "--pipeline" flag and add the "--output" flag to name the output folder. 

```
nextflow run main_AMR++.nf --pipeline eval_qc --output AMR++_results 
```


## Explore QC results
```
ls AMR++_results/
ls AMR++_results/QC_analysis/
ls AMR++_results/QC_analysis/FastQC/
ls AMR++_results/QC_analysis/FastQC/*
```

Now, let's look at the multiQC results which aggregates all these files into a single report.

# Run trim_qc pipeline

Input files for trim_qc:
* sample reads (`--reads` = ${baseDir}/data/raw/*_R{1,2}.fastq.gz)

Relevant default parameters:
* `--adapters` = "${baseDir}/data/adapters/nextera.fa"
* `--leading` = 3
* `--trailing` = 3
* `--slidingwindow` = "4:15"
* `--minlen = 36`

Output files:
* Trimmed reads

## trim_qc command

``` 
nextflow run main_AMR++.nf --pipeline trim_qc --output AMR++_results
```


### View trimmomatic.stats file on GUI
```
ls AMR++_results/QC_trimming/

ls AMR++_results/Results/Stats/
```

### Now, we'll talk about removing contaminant host DNA


---------------------

# Run host removal pipeline


Input files for `--pipeline rm_host` (need to change):
* QC trimmed reads (by default `--reads` = ${baseDir}/data/raw/*_R{1,2}.fastq.gz)

Relevant default parameters:
* `--host` = "${baseDir}/data/host/chr21.fasta.gz"

Output files:
* Non-host reads


## rm_host command

``` 
nextflow run main_AMR++.nf --pipeline rm_host --output AMR++_results --reads "AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz"
```

### Evaluate rm_host results
```
ls AMR++_results/HostRemoval/NonHostFastq/
ls AMR++_results/Results/Stats/
```

--------------------


# Run resistome pipeline

Input files for `--pipeline resistome`:
* `--reads`
* `--amr` = "${baseDir}/data/amr/megares_database_v3.00.fasta"
* `--annotation` = "${baseDir}/data/amr/megares_annotations_v3.00.csv"

Relevant default parameters:
* `--threshold` = 80
* `--min` = 5
* `--max` = 100
* `--skip` = 5
* `--samples` = 1

Optional parameters:
* `--snp` = "N"
* `--deduped` = "N"

Output files:
* Resistome analytic count matrix 
* Rarefaction analysis and plot

## resistome commmand

```
nextflow run main_AMR++.nf --pipeline resistome --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz"
```

## Evalute resistome results
```
ls AMR++_results/ResistomeAnalysis/
ls AMR++_results/ResistomeAnalysis/Rarefaction/Figures/
ls AMR++_results/Results/
```


----
# Run kraken pipeline

Input files for `--pipeline rm_host`:
* QC trimmed, nonhost sample reads (by default `--reads` = ${baseDir}/data/raw/*_R{1,2}.fastq.gz) 
* `--kraken_db` = null

Relevant default parameters:
* `--kraken_confidence` = 0.0

Output files:
* Kraken raw output
* Kraken classification reports
* Aggregated kraken results for all samples

## kraken command

``` 
nextflow run main_AMR++.nf --pipeline kraken --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" --kraken_db "/home/training/kraken2_DB/minikraken_8GB_20200312/"
```

## Inspect kraken results

```
ls 
```


----
# Effect of changing important parameters

Finally, let's explore how changing some of these parameters could affect our results. The two main examples are the `--threshold` flag which controls the gene fraction threshold in the resistome analysis, and the `--kraken_confidence` flag which controls how kraken classifies reads. 

```
nextflow run main_AMR++.nf --pipeline kraken --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" --kraken_db "/home/training/kraken2_DB/minikraken_8GB_20200312/" --kraken_confidence 1
```
