#!/bin/bash -ue
mkdir S2_test_fastqc_logs
fastqc -o S2_test_fastqc_logs -f fastq -q S2_test_R1.fastq.gz S2_test_R2.fastq.gz
