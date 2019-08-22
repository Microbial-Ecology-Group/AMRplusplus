Overview
--------

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AmrPlusPlus can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth. AmrPlusPlus can perform this analysis as well.

With AmrPlusPlus, you will obtain count files for each sample that can be combined into a count matrix and analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.


## Quick setup and test
First, make sure that Singularity is installed and in your $PATH variable. 
Visit this website for further information:
https://singularity.lbl.gov/docs-installation

### Download AmrPlusPlus v2
> git clone https://github.com/meglab-metagenomics/amrplusplus_v2.git
> cd amrplusplus_v2

### Install nextflow
> curl -s https://get.nextflow.io | bash

### Download minikraken database
> sh download_minikraken.sh

### Run test
> nextflow run main_amr_plus_plus_v2.nf -profile singularity


More Information
----------------

- [Software Requirements](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/requirements.md)
- [Installation](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/installation.md)
- [Usage](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/usage.md)
- [Configuration](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/configuration.md)
- [Output](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/output.md)
- [Dependencies](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/dependencies.md)
- [Contact](https://github.com/EnriqueDoster/bioinformatic-nextflow-pipelines/blob/master/docs/contact.md)
