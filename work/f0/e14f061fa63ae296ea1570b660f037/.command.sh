#!/bin/bash -ue
mkdir S3_test_fastqc_logs
fastqc -o S3_test_fastqc_logs -f fastq -q S3_test_R1.fastq.gz S3_test_R2.fastq.gz
