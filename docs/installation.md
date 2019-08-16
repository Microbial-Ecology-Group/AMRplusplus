Installation
------------

This section will help you get started with running the AmrPlusPlus pipeline with Nextflow and Docker. This tutorial assumes you will be running the pipeline from a POSIX compatible system such as Linux, Solaris, or OS X.

Setup
-----

We will go over a typical pipeline setup scenario in which you connect to a remote server, install Nextflow, and download the pipeline source code.

```bash
# username and host address
$ ssh [USER]@[HOST]

# install Nextflow
$ curl -s https://get.nextflow.io | bash

# give write permissions to user
$ chmod u+x nextflow

# move nextflow executable to a folder in your PATH environment variable
$ mv nextflow $HOME/bin

# create a test directory and change into it
$ mkdir amr_test && cd amr_test

# clone pipeline source code
$ git clone https://github.com/cdeanj/amrplusplus .
```

Run a Simple Test
-----------------

We will run a small sample dataset that comes with the pipeline source code. As such, we will not be specifying any input paths as they have already been included. During the program's execution, the required tool dependencies will be installed with Docker if they haven't been already. As there are many tool dependencies, this could take some time depending on your connection speed.

```bash
# command to run the amrplusplus pipeline
$ nextflow run main.nf -profile docker --threads 4 --output test

# change directories to view pipeline outputs
$ cd test/
```


