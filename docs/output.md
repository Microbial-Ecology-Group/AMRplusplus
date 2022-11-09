Output
------

All intermediate outputs produced from each module of this pipeline are provided as flat files that can be viewed in a text editor. These files are copied from the root **work/** directory created by Nextflow, so if disk space is a concern, this directory should be deleted as it can get quite large.

Directory Structure
-------------------

The output directories created by the pipeline are named after the module that produced them. Each file output is prefixed with the sample name and suffixed with a short product description. 

Files without sample prefixes are a result of aggregation. For example, the files **host.removal.stats** and **trimmomatic.stats** provide count matrices for the number of reads discarded as a result of host-dna removal and number of trimmed reads for each sample. 

Below is an example of all of the results created with a test run of AMR++.

```bash
test_results
├── Alignment
│   ├── BWA_Index
│   │   ├── chr21.fasta.gz.amb
│   │   ├── chr21.fasta.gz.ann
│   │   ├── chr21.fasta.gz.bwt
│   │   ├── chr21.fasta.gz.pac
│   │   ├── chr21.fasta.gz.sa
│   │   ├── megares_database_v3.00.fasta.amb
│   │   ├── megares_database_v3.00.fasta.ann
│   │   ├── megares_database_v3.00.fasta.bwt
│   │   ├── megares_database_v3.00.fasta.pac
│   │   └── megares_database_v3.00.fasta.sa
│   └── SAM_files
│       ├── Deduped
│       │   ├── S1_test.alignment.dedup.sam
│       │   ├── S2_test.alignment.dedup.sam
│       │   └── S3_test.alignment.dedup.sam
│       └── Standard
│           ├── S1_test.alignment.sam
│           ├── S2_test.alignment.sam
│           └── S3_test.alignment.sam
├── HostRemoval
│   └── NonHostFastq
│       ├── S1_test.non.host.R1.fastq.gz
│       ├── S1_test.non.host.R2.fastq.gz
│       ├── S2_test.non.host.R1.fastq.gz
│       ├── S2_test.non.host.R2.fastq.gz
│       ├── S3_test.non.host.R1.fastq.gz
│       └── S3_test.non.host.R2.fastq.gz
├── QC_analysis
│   ├── FastQC
│   │   ├── S1_test_fastqc_logs
│   │   │   ├── S1_test_R1_fastqc.html
│   │   │   ├── S1_test_R1_fastqc.zip
│   │   │   ├── S1_test_R2_fastqc.html
│   │   │   └── S1_test_R2_fastqc.zip
│   │   ├── S2_test_fastqc_logs
│   │   │   ├── S2_test_R1_fastqc.html
│   │   │   ├── S2_test_R1_fastqc.zip
│   │   │   ├── S2_test_R2_fastqc.html
│   │   │   └── S2_test_R2_fastqc.zip
│   │   └── S3_test_fastqc_logs
│   │       ├── S3_test_R1_fastqc.html
│   │       ├── S3_test_R1_fastqc.zip
│   │       ├── S3_test_R2_fastqc.html
│   │       └── S3_test_R2_fastqc.zip
│   └── MultiQC_stats
│       └── multiqc_report.html
├── QC_trimming
│   ├── Paired
│   │   ├── S1_test.1P.fastq.gz
│   │   ├── S1_test.2P.fastq.gz
│   │   ├── S2_test.1P.fastq.gz
│   │   ├── S2_test.2P.fastq.gz
│   │   ├── S3_test.1P.fastq.gz
│   │   └── S3_test.2P.fastq.gz
│   └── Unpaired
│       ├── S1_test.1U.fastq.gz
│       ├── S1_test.2U.fastq.gz
│       ├── S2_test.1U.fastq.gz
│       ├── S2_test.2U.fastq.gz
│       ├── S3_test.1U.fastq.gz
│       └── S3_test.2U.fastq.gz
├── ResistomeAnalysis
│   ├── Rarefaction
│   │   └── Counts
│   │       ├── S1_test.class.tsv
│   │       ├── S1_test.gene.tsv
│   │       ├── S1_test.group.tsv
│   │       ├── S1_test.mech.tsv
│   │       ├── S1_test.type.tsv
│   │       ├── S2_test.class.tsv
│   │       ├── S2_test.gene.tsv
│   │       ├── S2_test.group.tsv
│   │       ├── S2_test.mech.tsv
│   │       ├── S2_test.type.tsv
│   │       ├── S3_test.class.tsv
│   │       ├── S3_test.gene.tsv
│   │       ├── S3_test.group.tsv
│   │       ├── S3_test.mech.tsv
│   │       └── S3_test.type.tsv
│   └── ResistomeCounts
│       ├── S1_test.AMR.class.tsv
│       ├── S1_test.AMR.gene.tsv
│       ├── S1_test.AMR.group.tsv
│       ├── S1_test.AMR.mechanism.tsv
│       ├── S1_test.AMR.type.tsv
│       ├── S1_test.dedup_AMR.class.tsv
│       ├── S1_test.dedup_AMR.gene.tsv
│       ├── S1_test.dedup_AMR.group.tsv
│       ├── S1_test.dedup_AMR.mechanism.tsv
│       ├── S1_test.dedup_AMR.type.tsv
│       ├── S2_test.AMR.class.tsv
│       ├── S2_test.AMR.gene.tsv
│       ├── S2_test.AMR.group.tsv
│       ├── S2_test.AMR.mechanism.tsv
│       ├── S2_test.AMR.type.tsv
│       ├── S2_test.dedup_AMR.class.tsv
│       ├── S2_test.dedup_AMR.gene.tsv
│       ├── S2_test.dedup_AMR.group.tsv
│       ├── S2_test.dedup_AMR.mechanism.tsv
│       ├── S2_test.dedup_AMR.type.tsv
│       ├── S3_test.AMR.class.tsv
│       ├── S3_test.AMR.gene.tsv
│       ├── S3_test.AMR.group.tsv
│       ├── S3_test.AMR.mechanism.tsv
│       ├── S3_test.AMR.type.tsv
│       ├── S3_test.dedup_AMR.class.tsv
│       ├── S3_test.dedup_AMR.gene.tsv
│       ├── S3_test.dedup_AMR.group.tsv
│       ├── S3_test.dedup_AMR.mechanism.tsv
│       └── S3_test.dedup_AMR.type.tsv
└── Results
    ├── AMR_analytic_matrix.csv
    ├── SNPconfirmed_AMR_analytic_matrix.csv
    ├── SNPconfirmed_dedup_AMR_analytic_matrix.csv
    ├── Stats
    │   ├── host.removal.stats
    │   └── trimmomatic.stats
    └── dedup_AMR_analytic_matrix.csv



```
