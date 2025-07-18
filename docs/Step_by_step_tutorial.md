# Tutorial to run AMR++ in steps

The tutorial below takes you through the process of installing AMR++ and analyzing your data in steps, instead of running the entire pipeline. This can be useful if you wish to only perform certain components or your dataset size is extremely large and you need to reduce the size of the temporary files that are created in the working directory. 

It's important to note that nextflow populates the "work" directory with temporary files that allow for nextflow to "-resume" a failed run and pick up where it left off. Therefore, running AMR++ in steps allows for the removal of the "work" directory between every run. 

# Table of Contents
- [Tutorial to run AMR++ in steps](#tutorial-to-run-amr-in-steps)
- [Table of Contents](#table-of-contents)
- [Load conda environment](#load-conda-environment)
    - [Check conda installation](#check-conda-installation)
      - [Using server modules](#using-server-modules)
      - [Using miniconda](#using-miniconda)
    - [Install AMR++ environment](#install-amr-environment)
- [Run first Demo](#run-first-demo)
  - [How did AMR++ know what parameters to run:](#how-did-amr-know-what-parameters-to-run)
- [AMR++ components](#amr-components)
- [Run the eval\_qc pipeline](#run-the-eval_qc-pipeline)
  - [Run first pipeline, eval\_qc](#run-first-pipeline-eval_qc)
  - [Explore QC results](#explore-qc-results)
- [Run trim\_qc pipeline](#run-trim_qc-pipeline)
  - [trim\_qc command](#trim_qc-command)
    - [Download trimmomatic.stats file and open on excel](#download-trimmomaticstats-file-and-open-on-excel)
    - [Now, we'll talk about removing contaminant host DNA](#now-well-talk-about-removing-contaminant-host-dna)
- [Run host removal pipeline](#run-host-removal-pipeline)
  - [rm\_host command](#rm_host-command)
    - [Evaluate rm\_host results](#evaluate-rm_host-results)
- [Run resistome pipeline](#run-resistome-pipeline)
  - [resistome commmand](#resistome-commmand)
  - [Evalute resistome results](#evalute-resistome-results)
- [Run kraken pipeline](#run-kraken-pipeline)
  - [kraken command](#kraken-command)
  - [Inspect kraken results](#inspect-kraken-results)
- [Effect of changing important parameters](#effect-of-changing-important-parameters)
  - [Resistome](#resistome)
  - [Kraken](#kraken)


# Load conda environment
First, we have to have Anaconda or "conda" available for use. Then, we can install or load the conda environment, AMR++_env, which contains all the tools we need.

### Check conda installation
Try just running conda to check if it works for you:
```bash
conda --h
```

#### Using server modules

If it says `bash: conda: command not found...`, you have a few options to setup conda to work with your environment. In our university computing environment, modules are added that allow for easy access to other tools. We can run the following command:

```bash
module load Anaconda3/2024.02-1
```

This should load conda for you, but if it's the first time you ever do this, you'll have to run these two commands:
```bash
conda init

source ~/.bashrc
```

#### Using miniconda

If you don't have acess to conda, we recommend downloading [miniconda](https://docs.anaconda.com/miniconda/install/). This works the same as "conda" but comes in a smaller package for easier installation. 

Here are the basic commands for installation:
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
After installing, close and reopen your terminal application or refresh it by running the following command:

```bash
source ~/miniconda3/bin/activate
conda init --all
```

### Install AMR++ environment

Now, you should have conda installed and can install the AMR++ environment.

We'll need to download the code for [AMR++ from github](https://github.com/Microbial-Ecology-Group/AMRplusplus/tree/master) using the "git" command. 

```bash
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git
```

Let's navigate into the AMRplusplus directory and check the contents

```bash
cd AMRplusplus/

ls
```
We have a recipe to create the conda environment in the "envs" directory. 

```bash
conda env create -f envs/AMR++_env.yaml
# This can take 5-10 mins (or more) depending on your internet speed, computing resources, etc. 

# Once it's completed, activate the environment
conda activate AMR++_env

# You now have access to all the AMR++ software dependencies (locally)
samtools --help
```


# Run first Demo

We'll run the demonstration which will create output in the "test_results" directory. 

We can run the demo with this simple command:
```bash
nextflow run main_AMR++.nf
```

Note, we don't have to specify what `--profile` to use because it defaults to using the "local" profile and just looking in your regular $PATH, whis is now modified by the conda environment, AMR++_env.

## How did AMR++ know what parameters to run:

AMR++ has default parameters that are listed in the `params.config` file. We can then either change the values in the file directly, or add those flags to the command you're using directly. Also, note that there are variables with a single dash "-" and others with two dashes "--". The single dashes are internal to nextflow and include parameters/flags/variables like "-profile" and "-resume". The profile flag will not change if you're using conda as instructed above, and the "-resume" flag can be added to commands when the initial run failed for some reason and you want to try picking up where it left off. Otherwise, the majority of important variables will be denoted by two dashes "--". We'll go over all the parameters used by AMR++ which we'll either change in the command or by modifying the `params.config` file.  

Nextflow prioritizes:
1. whatever flags you include in the command
2. The default paramaters in `params.config`
3. Other variable calls within AMR++

We'll practice changing various variables below. 

First, let's look at the defaults:
```bash
less params.config 
```

For example, you can see which `--reads` are being analyzed by default and their location. We'll talk about each set of relevant variables as we go along.  

# AMR++ components

If the demo above completed succesfully, we are now ready run each component of the AMR++ pipeline by changing the `--pipeline` flag. Here is the list of pipelines included in AMR++. 

    Available pipelines:
        - demo: Run a demonstration of AMR++
        - standard_AMR: Run the standard AMR++ pipeline
        - fast_AMR: Run the fast AMR++ pipeline without host removal.
        - standard_AMR_wKraken: Run the standard AMR++ pipeline with Kraken
    Available pipeline subworkflows:
        - eval_qc: Run FastQC analysis
        - trim_qc: Run trimming and quality control
        - rm_host: Remove host reads
        - resistome: Perform resistome analysis
        - align: Perform alignment to MEGARes database
        - kraken: Perform Kraken analysis
        - qiime2: Perform QIIME 2 analysis
        - bam_resistome: Perform resistome analysis on BAM files

    To run a specific pipeline/subworkflow, use the "--pipeline" option followed by the pipeline name:
        nextflow run main_AMR++.nf --pipeline <pipeline_name> [other_options]

# Run the eval_qc pipeline

Input files for `--pipeline eval_qc`:
* sample reads (`--reads`) (by default `--reads` = ${baseDir}/data/raw/*_R{1,2}.fastq.gz)

Output files:
* FastQC results for each sample
* MultiQC report with aggregated results

## Run first pipeline, eval_qc

For this subworkflow, we'll just specify the "--pipeline" flag and add the "--output" flag to name the output folder. 

```bash
nextflow run main_AMR++.nf --pipeline eval_qc --output AMR++_results 
```


## Explore QC results
```bash
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

```bash
nextflow run main_AMR++.nf --pipeline trim_qc --output AMR++_results
```


### Download trimmomatic.stats file and open on excel
```bash
ls AMR++_results/QC_trimming/

ls AMR++_results/Results/Stats/
```

Don't forget to delete the "work" directory.

```bash
rm -r work/
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

```bash
nextflow run main_AMR++.nf --pipeline rm_host --output AMR++_results --reads "AMR++_results/QC_trimming/Paired/*{1,2}P.fastq.gz"
```

### Evaluate rm_host results
```bash
ls AMR++_results/HostRemoval/NonHostFastq/
ls AMR++_results/Results/Stats/
```

Don't forget to delete the "work" directory.

```bash
rm -r work/
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

```bash
nextflow run main_AMR++.nf --pipeline resistome --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz"
```

## Evalute resistome results
```bash
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

```bash 
nextflow run main_AMR++.nf --pipeline kraken --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" --kraken_db "/scratch/group/vero_research/databases/kraken2/PlusPF_8/"
```

## Inspect kraken results

```bash
ls 
```
Don't forget to delete the "work" directory.

```bash
rm -r work/
```

----
# Effect of changing important parameters

Finally, let's explore how changing some of these parameters could affect our results. The two main examples are the `--threshold` flag which controls the gene fraction threshold in the resistome analysis, and the `--kraken_confidence` flag which controls how kraken classifies reads. 


## Resistome 

We can rename the original reistome output AMR_analytic_matrix.csv to AMR_analytic_matrix_thresh80.csv to store the default results. Then, we can run the analysis again using a lower threshold. 

```bash
mv AMR++_results/Results/AMR_analytic_matrix.csv AMR++_results/Results/AMR_analytic_matrix_thresh80.csv

nextflow run main_AMR++.nf --pipeline resistome --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" --threshold 30

wc -l AMR++_results/Results/AM*
```

Notice the increase in taxa identified with the lower threshold.


## Kraken
Here, we can change the kraken_confidence score and run the kraken classification again. 

```bash
nextflow run main_AMR++.nf --pipeline kraken --output AMR++_results --reads "AMR++_results/HostRemoval/NonHostFastq/*R{1,2}.fastq.gz" --kraken_db "/scratch/group/vero_research/databases/kraken2/PlusPF_8/" --kraken_confidence 1

head AMR++_results/Results/unclassifieds_kraken_analytic_matrix.conf_*
```

Notice the major change in results, with a confidence of 1 leading to 100% unclassified reads.
