Installation
------------

This section will help you get started with running the AmrPlusPlus pipeline with Nextflow and Docker. This tutorial assumes you will be running the pipeline from a POSIX compatible system such as Linux, Solaris, or OS X.

Setup
-----

We will go over a typical pipeline setup scenario in which you connect to a remote server, install Nextflow, and download the pipeline source code. For the easist use of AmrPlusPlus, make sure that Singularity is installed and in your $PATH variable. 
Visit this website for further information:
https://singularity.lbl.gov/docs-installation

If Singularity cannot be installed, configure the "config/local.config" file to specify the absolute PATH to each required bioinformatic tool. Then, change the flag after "-profile" to "local" when running the pipeline.

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
$ chmod u+x nextflow

# move nextflow executable to a folder in your PATH environment variable
$ mv nextflow $HOME/bin

# create a test directory and change into it
$ mkdir amr_test && cd amr_test

# clone pipeline source code
$ git clone https://github.com/meglab-metagenomics/amrplusplus_v2.git .
```

Run a Simple Test
-----------------

We will run a small sample dataset that comes with the pipeline source code. As such, we will not be specifying any input paths as they have already been included. During the program's execution, the required tool dependencies will be accessed using a Singularity container. As there are many tool dependencies, downloading the container could take some time depending on your connection speed.

```bash
# navigate into AmrPlusPlus repository
$ cd amrplusplus_v2/

# command to run the amrplusplus pipeline
$ nextflow run main_amr_plus_plus_v2.nf -profile singularity --output test

# change directories to view pipeline outputs
$ cd test/
```


