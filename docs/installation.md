# Installation overview
This section will help you get started with running the AMR++ pipeline. This tutorial assumes you will be running the pipeline from a POSIX compatible system such as Linux, Solaris, or OS X.

There are four main ways that you can install and run AMR++, based on what is easiest for your computing cluster and whether conda, singularity, or docker is already installed. 

Based on our experience, we recommend following the instructions to [install miniconda without "sudo" permissions](#installing-miniconda-without-sudo-permissions).

Usually, you'll have to install nextflow, unless it's available to be loaded as module in your HPC.
* [Install Nextflow](#installing-nextflow)

To make all of the bioinformatic tool dependencies available for use with AMR++, we have a few options:
- [Installation overview](#installation-overview)
  - [Installing nextflow](#installing-nextflow)
  - [Run AMR++ with anaconda](#run-amr-with-anaconda)
    - [Installing miniconda without "sudo" permissions.](#installing-miniconda-without-sudo-permissions)
  - [Run AMR++ using Singularity](#run-amr-using-singularity)
  - [Run AMR++ using Apptainer](#Run-amr-using-apptainer)
  - [Run AMR++ using Docker](#run-amr-using-docker)
  - [Local installation of tools](#local-installation-of-tools)

## Installing nextflow

```bash
# username and host address
$ ssh [USER]@[HOST]

# Check if you have nextflow installed,
$ nextflow -h

# If not available, install Nextflow
$ curl -s https://get.nextflow.io | bash
# If you do not have curl installed, try wget
# $ wget -qO- https://get.nextflow.io | bash

# give write permissions to user
chmod u+x nextflow

# move nextflow executable to a folder in your $PATH environment variable. For example:
mv nextflow $HOME/bin
```

## Run AMR++ with anaconda
Requirements:
* Nextflow
* Anaconda

If anaconda is already installed and nextflow is working, we'll just need to download the AMR++ github repository.

```bash
# Download AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

# Navigate into direcotry
cd AMRplusplus

# Install mamba for faster installation
conda install mamba -n base -c conda-forge

# Test AMR++ by specifying the "conda" profile. 
nextflow run main_AMR++.nf -profile conda

# If your computing cluster uses the slurm scheduler, modify the script "run_AMR++_slurm.sh" to
# accurately request computing resources. Then, run it using:
sbatch run_AMR++_slurm.sh
```



### Installing miniconda without "sudo" permissions.

We will go over a typical pipeline setup scenario in which you connect to a remote server, install miniconda (or use a local installation of anaconda), and download the pipeline source code. In cases where the Anaconda installation on your computing cluster is not updated or you are experiencing errors while installing packages, we recommend miniconda. Use [this site](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) for further information on installing miniconda on a user-writable directory. 

```bash
# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh
# Run installation, follow default options. However, depending on your computing cluster when you are prompted, you should consider changing the 
# location of the installation to somewhere that is not your home directory which can have storage limits.
bash Miniconda3-py310_23.11.0-2-Linux-x86_64.sh

# Download AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

# Navigate into direcotry
cd AMRplusplus

# Test AMR++ by specifying the "conda" profile. If your computing cluster uses the slurm scheduler, use the "conda_slurm" profile.
nextflow run main_AMR++.nf -profile conda
```

Optionally, you can create the conda environment prior to running AMR++.

```bash
# If you haven't created the conda AMR++ environment as instructed above, go ahead and do that here:
# If you have mamba installed, you can swap "conda" with "mamba" for faster installation.
conda env create -f envs/AMR++_env.yaml
# Now activate the conda environment
conda activate activate AMR++_env

# Now, the tools will be available "locally" so we must run AMR++ using the "local" profile.
nextflow run main_AMR++.nf -profile local

```

## Run AMR++ using Singularity

Requirements:
* Nextflow
* Singularity

Often, HPCs will have singularity installed this will allow AMR++ to download and use a singularity container with all of the pre-installed software requirements.

```bash

# Download AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

# Navigate into direcotry
cd AMRplusplus

# Run command with singularity profile
nextflow run main_AMR++.nf -profile singularity

# Alternatively, you can pull the singularity container first like this:
singularity pull docker://enriquedoster/amrplusplus:latest

# Then, specify the path to the singularity image.
nextflow run main_AMR++.nf -profile local -with-singularity amrplusplus_latest.sif

```

## Run AMR++ using Apptainer

Requirements:
* Nextflow
* Apptainer

Apptainer is the opensourc fork of singularity that is a part of the linux project. Sometimes HPCs might have apptainer intstalled rather than singularity, this will allow AMR++ to download and use a singularity/apptainer container with all of the pre-installed software requirements.

```bash

# Download AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

# Navigate into direcotry
cd AMRplusplus

# Run command with singularity profile
nextflow run main_AMR++.nf -profile apptainer

# Alternatively, you can pull the singularity container first like this:
singularity pull docker://enriquedoster/amrplusplus:latest

# Then, specify the path to the singularity image.
nextflow run main_AMR++.nf -profile local -with-apptainer amrplusplus_latest.sif

```

## Run AMR++ using Docker

Requirements:
* Nextflow
* Docker

Like Singularity, Docker is another tool management option that is often available on HPCs and this can be used by AMR++ to download a docker container with all of the pre-installed software requirements.

```bash

# Download AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

# Navigate into directory
cd AMRplusplus

# Run command with docker profile
nextflow run main_AMR++.nf -profile docker

```



## Local installation of tools

Requirements:
* All [software dependencies](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/dependencies.md)
* Nextflow


If you run into issues with using Miniconda and Singularity, or perhaps your computing cluster has all of the required tools, configure the "config/local.config" file to specify the absolute PATH to each required bioinformatic tool. Or if you can add all of the relevant paths to your $PATH environment variable, which in some cases means loading the correct modules, then you can just leave the name of each tool without the path. Finally, change the flag after "-profile" to "local" when running the pipeline.

```bash
nextflow run main_AMR++.nf -profile local
```
