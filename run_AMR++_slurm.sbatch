#!/bin/bash
#SBATCH -J AMR++ -o AMR++_log.out -t 24:00:00 --mem=5G --nodes=1 --ntasks=1 --cpus-per-task=1

# Remember to change the "maxForks" parameter in the config file that you are using. This corresponds with "--ntasks" 
# to control how many jobs are submitted at once. The "--threads" argument should also match the --cpus-per-task to 
# fully utlize the available computing resources.

# This script works on TAMU's Grace HPRC, but you need to follow the instructions on the Github to get the right conda 
# environment installed on your computing environment
conda activate AMR++_env  # Explore the installation instructions on github to see how to install this environment

# If you are submitting this script using sbatch to submit to the SLURM scheduler, we use the "-profile local_slurm". 
# This profile will submit individual jobs to SLURM with the default resources detailed in "config/local_slurm.config" file.
nextflow run main_AMR++.nf -profile local_slurm --threads 8 


