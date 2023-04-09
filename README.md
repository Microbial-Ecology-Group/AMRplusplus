Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)


# AMR++ bioinformatic pipeline
(https://megares.meglab.org/)

AMR++ is a bioinformatic pipeline meant to aid in the analysis of raw sequencing reads to characterize the profile of antimicrobial resistance genes, or resistome. AMR++ was developed to work in conjuction with the the MEGARes database which contains sequence data for approximately 9,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AMR++ can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth.

With AMR++, you will obtain alignment count files for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

More Information
----------------

- [Installation](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/installation.md)
- [Usage](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/usage.md)
- [Configuration](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/configuration.md)
- [Output](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/output.md)
- [Dependencies](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/dependencies.md)
- [Software Requirements](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/requirements.md)
- [FAQs](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/FAQs.md)
- [Details on AMR++ updates](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/update_details.md)
- [Contact](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/contact.md)



## AMR++ demonstration

If anaconda is already installed and nextflow is working, we'll just need to download the AMR++ github repository. Please review the [installation document](docs/installation.md) for alternative methods to install AMR++ in your computing environment.

```bash
# Install mamba for faster installation
conda install mamba -n base -c conda-forge
```

Clone the AMR++ repository.

```bash
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git
```

Navigate into the AMR++ repository and run the test command.
```bash
cd AMRplusplus

# Run command to perform the demonstration pipeline using the conda profile.
nextflow run main_AMR++.nf -profile conda

# The first time this can take 5-10 mins (or more) depending on your internet speed because it is installing a conda environment. Subsequent runs will skip this step automatically.
```
Now, you can check out the results in the newly created "test_results" directory.

# Using AMR++ to analyze your data

AMR++ is customizable to suit your computing needs and analyze your data. Primarily, the ```-profile``` paramater allows you to choose between running AMR++ using a singularity container, docker container, anaconda packages, or a local installation of your software. 
All parameters used to control how AMR++ analyzes your data can also be changed as needed in a variety of ways. For full information, review this [configuration document.](docs/configuration.md)


Below is a brief example, the default parameters were run using this command:

```nextflow run main_AMR++.nf -profile conda```

To change the reads that were analyzed, you should specify the ```--reads`` parameters. Here, we can use regular expressions to point to your samples in a different directory.
```bash
nextflow run main_AMR++.nf -profile conda --reads "path/to/your/reads/*_R{1,2}.fastq.gz" 
```

# Optional flags

## SNP verification

AMR++ now works in conjuction with a [custom SNP verification software](https://github.com/Isabella136/AmrPlusPlus_SNP) to evaluate alignments to gene accessions requiring SNP confirmation to confer resistance. To include this workflow, include the ```--snp Y``` flag in your command like this:

```bash
nextflow run main_AMR++.nf -profile conda --snp Y
```
This will create with the standard count table (AMR_analytic_matrix.csv) in addition to a count matrix with SNP confirmed counts (SNPconfirmed_AMR_analytic_matrix.csv).

## Deduplicated counts

Another option is to include results for deduplicated counts by using the ```--deduped Y``` flag in your command.

```bash
nextflow run main_AMR++.nf -profile conda --snp Y --deduped Y
```

With this flag, AMR++ will extract the deduplicated alignments to MEGARes also output a count matrix with deduplicated counts. Since also we included the ```--snp Y``` flag, we will end up with 4 total output count matrices.

# Choosing the right pipeline

AMR++ analyzes data by combining workflows that takes a set of sequencing reads through various bioinformatic software. We recommend our standard AMR++ pipeline as a comprehensive way to start from raw sequencing reads, QC assessment, host DNA removal, and resistome analysis with MEGARes. However, users might only want to replicate portions of the pipeline and have more control over their computing needs. Using the ```--pipeline``` parameter, users can now change how AMR++ runs.



## Pipeline workflows
*  omitting the ```--pipeline``` flag or using ```--pipeline demo```    
    * Simple demonstration on test data

* ```--pipeline standard_AMR```   
    * Steps: QC trimming > Host DNA removal > Resistome alignment > Resistome results

* ```--pipeline fast_AMR```
    * This workflow simply skips host removal to speed up analysis.
    * Steps: QC trimming > Resistome alignment > Resistome results

* ```--pipeline standard_AMR_wKraken```
    * This workflow adds microbiome analysis with kraken. It requires having a local kraken database. The minikraken_8GB_202003 will be downloaded automatically and requires ~8GB of space. Otherwise, you can specify the location to your own database with the flag, ```--kraken_db "/Path/to/KrakenDb/"```
    * Steps:
        * QC trimming > Host DNA removal > Resistome alignment > Resistome results 
        * Non-host reads > Microbiome analysis

## Pipeline subworkflows
* ```--pipeline eval_qc```  
    * Evaluate sample QC 
* ```--pipeline trim_qc```  
    * QC trimming using trimmomatic 
* ```--pipeline rm_host```  
    * Align reads to host DNA using bwa and remove contaminants 
* ```--pipeline resistome```  
    * Align reads to MEGARes using bwa, perform rarefaction and resistome analysis
* ```--pipeline kraken```  
    * Classify reads taxonomically using kraken.
* ```--pipeline bam_resistome```
    * This will run the resistome pipeline starting with bam files from a previous alignment to MEGARes.
    * Need to include ```--bam_files "Path/to/BAM/*.bam"``` in the command line.

## Example command
In the following example, we'll choose to run the standard AMR++ workflow, which includes QC trimming, host removal, and Resistome analysis. Since we included the ```--snp Y --deduped Y``` flags, we'll also get ouput for deduped counts and SNP confirmed counts.

Alternatively, you can modify all of these variables and more in the "params.config" file which will be loaded automatically. Just make sure to include the "-profile" and "--pipeline" flags. More information [in this document](docs/configuration.md)

```bash
# Remember to update the --reads flag to match your read location
nextflow run main_AMR++.nf -profile conda --pipeline standard_AMR --reads "path/to/your/reads/*_R{1,2}.fastq.gz" --snp Y --deduped Y
```
