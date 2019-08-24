Output
------

All intermediate outputs produced from each module of this pipeline are provided as flat files that can be viewed in a text editor. These files are copied from the root **work/** directory created by Nextflow, so if disk space is a concern, this directory should be deleted as it can get quite large.

Directory Structure
-------------------

The output directories created by the pipeline are named after the module that produced them. Each file output is prefixed with the sample name and suffixed with a short product description. 

Files without sample prefixes are a result of aggregation. For example, the files **host.removal.stats** and **trimmomatic.stats** provide count matrices for the number of reads discarded as a result of host-dna removal and number of trimmed reads for each sample. 

```bash
├── RunQC
│   ├── Paired
│   │   ├── SRR532663.1P.fastq
│   │   └── SRR532663.2P.fastq
│   ├── Stats
│   │   └── trimmomatic.stats
│   └── Unpaired
│       ├── SRR532663.1U.fastq
│       └── SRR532663.2U.fastq
├── BuildHostIndex
│   ├── chr21.fasta.amb
│   ├── chr21.fasta.ann
│   ├── chr21.fasta.bwt
│   ├── chr21.fasta.pac
│   └── chr21.fasta.sa
├── AlignReadsToHost
│   └── SRR532663.host.sam
├── NonHostReads
│   ├── SRR532663.non.host.R1.fastq
│   └── SRR532663.non.host.R2.fastq
├── RemoveHostDNA
│   ├── HostRemovalStats
│   │   └── host.removal.stats
│   └── NonHostBAM
│       └── SRR532663.host.sorted.removed.bam
├── AlignToAMR
│   └── SRR532663.amr.alignment.sam
├── RunResistome
│   └── SRR532663.gene.tsv
├── ResistomeResults
│   └── AMR_analytic_matrix.csv
├── RunRarefaction
│   ├── SRR532663.class.tsv
│   ├── SRR532663.gene.tsv
│   ├── SRR532663.group.tsv
│   └── SRR532663.mech.tsv
└── RunSNPFinder
│   └── SRR532663.tsv
└── RunFreebayes
│   └── SRR532663.results.vcf
├── RunKraken
│   └── SRR532663.kraken.report
│   └── SRR532663.kraken.filtered.report
├── KrakenResults
│   └── kraken_analytic_matrix.csv
├── FilteredKrakenResults
│   └── filtered_kraken_analytic_matrix.csv



```
