#!/bin/bash
#qsub -N=amrplusplus
#qsub -q small
#qsub -l procs-per-node 1
#qsub -l walltime 20:00:00
#qsub -l M user_email@gmail.com

module purge
module load jdk/1.8.0
module load singularity/2.5.2
module load ompi/gnu.mesabi

mpirun --pernode ./nextflow run main_AmrPlusPlus_v2.nf -resume -profile msi_pbs \
-w /work_dir --threads 15 \
--output output_results --host /PATH/TO/HOST/GENOME \
--reads "RAWREADS/*_R{1,2}.fastq.gz" -with-mpi
