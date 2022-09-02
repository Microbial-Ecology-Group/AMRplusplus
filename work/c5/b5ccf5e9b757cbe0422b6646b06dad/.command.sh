#!/bin/bash -ue
bwa mem -t 4 megares_modified_database_v2.00.fasta S3_test.1P.fastq.gz S3_test.2P.fastq.gz > S3_test.sam
samtools view -@ 4 -Sb > S3_test.bam
# LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
