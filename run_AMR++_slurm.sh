#!/bin/bash
#SBATCH -J AMR++ -o AMR++_log.out -t 24:00:00 -p amilan --nodes=1 --ntasks=4 --cpus-per-task=8

# Remember to change the "maxForks" parameter in the config file that you are using. This corresponds with "--ntasks" 
# to control how many jobs are submitted at once. The "--threads" argument should also match the --cpus-per-task to 
# fully utlize the available computing resources.

# This script works on TAMU's HPRC, but you need to follow the instructions on the Github to get the right conda 
# environment installed on your computing environment
module load Nextflow

conda activate AMR++_env  # Explore the installation instructions on github to see how to install this environment

nextflow run main_AMR++.nf -profile local_slurm --threads 8
