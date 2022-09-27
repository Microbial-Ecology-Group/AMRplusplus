Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)


# MEGARes and the AMR++ bioinformatic pipeline
(https://megares.meglab.org/)

The MEGARes database contains sequence data for approximately 9,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

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

Create the anaconda environment for AMR++. This will work for both the nextflow version and snakemake version.

```bash
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n AMR++ snakemake nextflow git make cxx-compiler singularity
mamba activate AMR++
```


Clone the AMR++ repository.

```bash
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git
```

Brief tutorial for nextflow pipeline test run
```bash
cd AMRplusplus

# Run command to perform the demonstration pipeline using the singularity profile
nextflow run main_AMR++.nf -profile singularity --pipeline demo
```


# Using AMR++ to analyze your data

AMR++ is customizable to suit your computing needs and analyze your data. Primarily, the ```-profile``` paramater allows you to choose between running AMR++ using a singularity container, docker container, anaconda packages, or a local installation of your software. 
All parameters used to control how AMR++ analyzes your data can also be changed as needed in a variety of ways. For full information, review this [configuration document.](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/configuration.md)


Below is a brief example, the default parameters were run using this command:

```nextflow run main_AMR++.nf -profile singularity --pipeline demo```

To change the reads that were analyzed, you should specify the ```--reads`` parameters. Here, we can use regular expressions to point to your samples in a different directory.

```nextflow run main_AMR++.nf -profile singularity --pipeline demo --reads "path/to/your/reads/*_R{1,2}.fastq.gz" ```


## Choosing a modified pipeline

AMR++ analyzes data by combining workflows that takes a set of sequencing reads through various bioinformatic software. We recommend our standard AMR++ pipeline as a comprehensive way to start from raw sequencing reads, QC assessment, host DNA removal, and resistome analysis with MEGARes. However, users might only want to replicate portions of the pipeline and have more control over their computing needs. Using the ```--pipeline``` parameter, users can change how AMR++ runs.

* ```--pipeline demo```    
    * Simple demonstration

* ```--pipeline standard_AMR```   
    * QC trimming > Host DNA removal > Resistome alignment > Resistome results

* ```--pipeline fast_AMR```
    * QC trimming > Resistome alignment > Resistome results

* ```--pipeline standard_AMR_wKraken```
    * QC trimming > Host DNA removal > Resistome alignment > Resistome results 
Non-host reads > Microbiome analysis

* pipeline fragments
    * ```--pipeline multiqc```  
        * Evaluate sample QC 
    * ```--pipeline trim```  
        * QC trimming using trimmomatic 
    * ```--pipeline rmhost```  
        * Align reads to host DNA using bwa and remove contaminants 
    * ```--pipeline resistome```  
        * Align reads to MEGARes using bwa, perform rarefaction and resistome analysis
    * ```--pipeline kraken```  
        * Classify reads taxanomically using kraken 

