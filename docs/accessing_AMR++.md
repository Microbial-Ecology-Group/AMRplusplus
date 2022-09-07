# Accessing AMR++



This section will help you get access to all the bioinformatic tools required for metagenomic analysis with AMR++.

## Anaconda environment
-----------------


## Miniconda

In cases where the Anaconda installation on your computing cluster is not updated or you are experiencing errors while installing packages, we recommend miniconda. Use [this site](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) for further information on installing miniconda on a user-writable directory. 


```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh

```


## Docker container
-----------------

Docker containers allow the packaging of multiple bioinformatic tools. While docker is a popular tool and likely to be supported by many computing clusters, please contact your system administrator for help with installing docker. Installation on a local computer is also an option and can be performed by following these instructions: https://docs.docker.com/get-docker/

We provide AMR++ with a docker container that is automatically accessed when running the AMR++ pipeline by using the flag, "-profile docker". 



