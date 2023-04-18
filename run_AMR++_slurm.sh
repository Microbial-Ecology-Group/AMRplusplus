#!/bin/bash
#SBATCH -J AMR++ -o AMR++_log.out -t 24:00:00 -p amilan --nodes=1 --ntasks=4 --cpus-per-task=8

# Remember to change the "queueSize" parameter in the config/*_slurm.config file that you are using. This corresponds with "--ntasks" 
# to control how many jobs are submitted at once.

# This script works on TAMU's HPRC, but you need to follow the instructions on the Github to get the right conda 
# environment installed on your computing environment
module load Nextflow

conda activate AMR++

nextflow run main_AMR++.nf -profile conda_slurm --threads 8
