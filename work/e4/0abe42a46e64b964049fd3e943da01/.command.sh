#!/bin/bash -ue
mkdir S1_test_fastqc_logs
fastqc -o S1_test_fastqc_logs -f fastq -q S1_test_R1.fastq.gz S1_test_R2.fastq.gz
