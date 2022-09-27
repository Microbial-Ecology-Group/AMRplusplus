Installation
------------

This section will help you get started with running the AMR++ pipeline. This tutorial assumes you will be running the pipeline from a POSIX compatible system such as Linux, Solaris, or OS X.

Setup
-----

We will go over a typical pipeline setup scenario in which you connect to a remote server, install miniconda (or use a local installation of anaconda), and download the pipeline source code. In cases where the Anaconda installation on your computing cluster is not updated or you are experiencing errors while installing packages, we recommend miniconda. Use [this site](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) for further information on installing miniconda on a user-writable directory. 

```bash
# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
# Run installation, follow default options. However, depending on your computing cluster when you are prompted, you should consider changing the 
# location of the installation to somewhere that is not your home directory which can have storage limits.
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh

# If you haven't created the mamba environment and AMR++ environment as instructed on the README file, go ahead and do that here:
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n AMR++ nextflow singularity
mamba activate AMR++



```

Run a Simple Test
-----------------

We will run a small sample dataset that comes with the pipeline source code. As such, we will not be specifying any input paths as they have already been included. During the program's execution, the required tool dependencies will be accessed using a Singularity container. As there are many tool dependencies, downloading the container could take some time depending on your connection speed.

```bash
# You are now ready to navigate into github repository
cd AMRplusplus

# Run test 
nextflow run main_AMR++.nf -profile singularity --pipeline demo

# change directories to view pipeline outputs
cd test_results/
```

## Installation without anaconda, using Singularity

Requirements:
* Nextflow
* Singularity

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

# move nextflow executable to a folder in your PATH environment variable
mv nextflow $HOME/bin

# Download AMR++ repository
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git

# Navigate into direcotry
cd AMRplusplus

# Run command with singularity profile
nextflow run main_AMR++.nf -profile singularity --pipeline demo

# Alternatively, you can pull the singularity container first like this:
singularity pull docker://enriquedoster/amrplusplus:latest

# 
nextflow run main_AMR++.nf -profile local --pipeline demo -with-singularity amrplusplus_latest.sif


```

## Local installation of tools

Requirements:
* All [software requirements](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/requirements.md)
* Nextflow


If Singularity cannot be installed, or perhaps your computing cluster has all of the required tools, configure the "config/local.config" file to specify the absolute PATH to each required bioinformatic tool. Or if you can add all of the relevant paths to your $PATH environment variable, which in some cases means loading the correct modules, then you can just leave the name of each tool without the path. Finally, change the flag after "-profile" to "local" when running the pipeline.

```nextflow run main_AMR++.nf -profile local --pipeline demo```