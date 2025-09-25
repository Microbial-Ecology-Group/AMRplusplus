#!/bin/bash
#SBATCH -J AMR++ -o AMR++_log.out -t 24:00:00 --mem=5G --nodes=1 --ntasks=1 --cpus-per-task=1

# Remember to change the "maxForks" parameter in the config file that you are using. This corresponds with "--ntasks" 
# to control how many jobs are submitted at once. The "--threads" argument should also match the --cpus-per-task to 
# fully utlize the available computing resources.

# This script works on TAMU's Grace HPRC, but you need to follow the instructions on the Github to get the right conda 
# environment installed on your computing environment
conda activate AMR++_env  # Explore the installation instructions on github to see how to install this environment

# Because we are using the local profile, software from the conda environment will be used. The flag "maxForks" is 
# in the local.config file and will spawn 4 processes by default, this corresponds with "--ntasks=4" in the sbatch script.
nextflow run main_AMR++.nf -profile local --threads 8 # This will use 8 threads, which corresponds with "--cpus-per-task=8". 

