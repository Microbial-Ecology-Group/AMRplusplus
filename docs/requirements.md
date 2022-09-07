Software Requirements
---------------------
To run AMR++, you will need the following tool installed on your server or local machine.

  - Anaconda or miniconda (Required)
    - Visit this website for further information: https://docs.anaconda.com/anaconda/install/
  
NOTE: If you choose not to install anaconda, you will need to download each of the required dependencies and add the executable paths to your .bashrc file to run the pipeline. A list of these dependencies can be found in the [Dependencies](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/dependencies.md) section of this document.

If anaconda or docker is not available to you, in addition to installing the listed bioinformatic tool dependencies, you'll need to install the following tools:

  - Java 7+ (Required)
  - Nextflow (Required)

Another option is to use singularity:

```bash
singularity pull docker://enriquedoster/amrplusplus:latest

nextflow run main_AMR++.nf -profile local --pipeline demo -with-singularity amrplusplus_latest.sif

```