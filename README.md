Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

## DRAFT ## DRAFT

## AMR++ v3 updates

Updates
- SNP confirmation software
- pipeline structure
 - switch to nextflow DSL2 and collection of "pipeline" options
 - snakemake approach for optimizing storage requirements
- Includes updates to MEGARes v3


Create the environment for AMR++. This will work for both the nextflow version and snakemake version.

```bash
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n AMR++ snakemake git nextflow
mamba activate AMR++
```


Clone the AMR++ repository.

```bash
mamba activate AMR++
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git
```

Brief tutorial for nextflow pipeline test run
```bash
cd AMRplusplus
mamba activate AMR++

nextflow run main_AMR++.nf -profile conda --pipeline demo
```


# Customizing the pipeline to analyze your data


## Changing the default variables

The pipeline comes with test datain the `data/` directory and uses default paramaters, found in the `params.config` file.
You can edit this file directly to change the parameters, or you can specify parameters on the command line using ``--<parameter name>```.

For example, the default parameters were run using this command:

```nextflow run main_AMR++.nf -profile conda --pipeline demo```

To change the reads that were analyzed, you should specify the ```--reads`` parameters. Here, we can use regular expressions to point to your samples in a different directory.

```nextflow run main_AMR++.nf -profile conda --pipeline demo --reads "path/to/your/reads/*_R{1,2}.fastq.gz" ```


## Profiles

```-profile conda``` 

```-profile local``` 

```-profile docker```


### Pipelines

```--pipeline demo```    Simple demonstration

```--pipeline standard_AMR```   QC trimming > Host DNA removal > Resistome alignment > Resistome results

```--pipeline fast_AMR```  QC trimming > Resistome alignment > Resistome results

```--pipeline standard_AMR_wKraken```   QC trimming > Host DNA removal > Resistome alignment > Resistome results 
Non-host reads > Microbiome analysis


pipeline fragments
```--pipeline multiqc```  Evaluate sample QC 




## Run SnakeMake Workflow

Brief tutorial for snakemake pipeline test run
```bash
cd AMRplusplus
mamba activate AMR++

snakemake --use-conda --cores <number of threads available>
```

# Microbial Ecology Group (MEG)
(https://megares.meglab.org/)

Our international multidisciplinary group of scientists and educators is addressing the issues of antimicrobial resistance (AMR) and microbial ecology in agriculture through research, outreach, and education. By characterizing risks related to AMR and microbial ecology, our center will identify agricultural production practices that are harmful and can be avoided, while also identifying and promoting production practices and interventions that are beneficial or do no harm to the ecosystem or public health. This will allow society to realize “sustainable intensification” of agriculture.

# MEGARes and the AMR++ bioinformatic pipeline
(http://megares.meglab.org/amrplusplus/latest/html/v2/)

The MEGARes database contains sequence data for approximately 9,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AmrPlusPlus can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth.

With AMR++, you will obtain alignment count files for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

More Information
----------------

- [Installation](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/installation.md)
- [Usage](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/usage.md)
- [Configuration](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/configuration.md)
- [Accessing AMR++](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/accessing_AMR++.md)
- [Output](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/output.md)
- [Dependencies](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/dependencies.md)
- [Software Requirements](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/requirements.md)
- [FAQs](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/FAQs.md)
- [Details on AMR++ updates](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/update_details.md)
- [Contact](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/contact.md)
