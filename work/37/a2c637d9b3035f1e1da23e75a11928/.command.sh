#!/bin/bash -ue
bwa mem -t 4 megares_modified_database_v2.00.fasta S1_test.1P.fastq.gz S1_test.2P.fastq.gz > S1_test.sam
samtools view -@ 4 -Sb > S1_test.bam
# LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
