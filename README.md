Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)


# AMR++ bioinformatic pipeline
(https://megares.meglab.org/)

AMR++ is a bioinformatic pipeline meant to aid in the analysis of raw sequencing reads to characterize the profile of antimicrobial resistance genes, or resistome. AMR++ was developed to work in conjuction with the the MEGARes database which contains sequence data for approximately 9,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AMR++ can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth.

With AMR++, you will obtain alignment count files for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

## Important changes to AMR++

[Detailed changes here.](docs/CHANGELOG.md)

Brief overview:
1. **Resistome counting now uses `alignment_analyzer.py` (AlignmentAnalyzer), replacing the compiled `resistomeanalyzer` binary.** The new analyzer reads the BAM directly and produces the gene-level count matrix, with several improvements to *how* alignments are turned into counts (items 2–4 below). The old `resistomeanalyzer` binary is no longer used for gene counting.
2. **Fragment-level counting — a major change to how hits are counted.** AMR++ now counts at the *fragment* level rather than counting every read as a separate hit.
   - Previously: every aligned read counted as one hit, so a paired-end fragment (R1 + R2) could be counted twice while a merged read counted once — inflating and skewing counts between library types.
   - Now: a read pair and a merged read each count as **one hit**, so paired-end and merged data are directly comparable. This is the default (`count_mode = "fragment"`).
   - When a fragment's two mates map to different-but-near-identical gene accessions in the same MEGARes **Group** (a common database-redundancy artifact), the hit is collapsed to one gene (`group_aware = "Y"`) rather than double-counted. Mates mapping to genuinely different Groups are counted once to each.
   - Query coverage is now "edge-aware" (`edge_aware_qcov = "Y"`): a read that aligns cleanly but runs off the end of a short gene isn't penalized for the overhanging portion. All three behaviors can be turned off (`count_mode = read_end`, `group_aware = "N"`, `edge_aware_qcov = "N"`).
3. **Only primary alignments are counted** for the resistome (secondary/supplementary alignments are excluded by default).
4. Changed default AMR gene fraction `--threshold` to `0`, and added query-coverage / gene-fraction filters to the analyzer (`--min_gene_fraction`, `--min_query_coverage`, both on a 0–1 proportion scale). We recommend running statistical analysis of count matrices after aggregating to the "Group" level to account for possible false-positive calls of individual gene accessions.
5. **Deduplication: alignment deduplication is no longer run by default — we now recommend read deduplication instead**, particularly when analyzing target-enriched data.
   - **Read deduplication** (`read_dedup = "Y"`, Clumpify, applied after QC trimming) more directly addresses PCR duplication by collapsing reads that are actual sequence duplicates.
   - **Alignment deduplication** (`deduped`) produced a larger overall reduction in read numbers, but that reduction included reads that were **not** true duplicates — reads could be collapsed based on alignment coordinates/length rather than being genuine PCR duplicates. Because of this over-collapsing, it is no longer the default recommendation.
   - Alignment deduplication is still available for comparison via `deduped = "Y"`; read deduplication is the recommended approach.
6. Added single-end and merged-read analysis.
7. Changed defaults to skip rarefaction analysis, but default to running SNP confirmation.
8. **Updated SNP confirmation tool (on a separate branch).** The SNP verifier now adjusts counts by the *proportion of reads that actually carry the confirmed SNP*: it computes (reads with the SNP / reads covering the SNP position) and multiplies by the gene's count, rather than substituting a raw resistant-read count. Genes that require SNP confirmation but have no entry in the SNP database are now set to zero (previously they passed through unconfirmed). The tool has also been updated to run cleanly under **Nextflow DSL2 syntax**.
9. **New coverage-threshold sweep workflow (`bam_coverage_sweep`).** Runs a set of pre-aligned BAMs across a grid of gene-fraction and query-coverage thresholds, then produces combined matrices and drop-off plots showing how the resistome (genes, classes, mechanisms, groups) shrinks as thresholds tighten — a way to choose sensible cutoffs and see each sample's sensitivity to them. See the [coverage sweep tutorial](docs/Coverage_sweep_tutorial.md).


[Additional analysis tips here.](docs/Analysis_recommendations.md)

More Information
----------------

- [Getting Started](docs/GettingStarted.md)
- [Installation](docs/installation.md)
- [Usage](docs/usage.md)
  - [Choosing the right pipeline](docs/choosing_pipeline.md)
  - [Analysis recommendations](docs/Analysis_recommendations.md)
  - [Paired-end analysis step-by-step](docs/Step_by_step_tutorial.md)
  - [Single-end analysis step-by-step](docs/SingleEnd_read_tutorial.md)
  - [Merged-read analysis step-by-step](docs/Merged_read_tutorial.md)
  - [Coverage threshold sweep](docs/Coverage_sweep_tutorial.md)
- [Configuration](docs/configuration.md)
  - [Tips for using SLURM](docs/Running_with_SLURM.md)
- [Output](docs/output.md)
- [Dependencies](docs/dependencies.md)
- [Troubleshooting and FAQs](FAQs.md)
- [Details on AMR++ updates](docs/CHANGELOG.md)
- [Contact](docs/contact.md)



## AMR++ demonstration

If anaconda is already installed, we'll just need to download the AMR++ github repository and create the AMR++ conda environment. Please review the [installation document](docs/installation.md) for alternative methods to install AMR++ in your computing environment.

```bash
# Confirm conda works
conda -h
```

Clone the AMR++ repository.

```bash
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus.git
```

Navigate into the AMR++ repository and run the test command.
```bash
cd AMRplusplus

# Now we can use the included recipe to make the AMR++ environment
conda env create -f envs/AMR++_env.yaml
# This can take 5-10 mins (or more) depending on your internet speed, computing resources, etc. 

# Once it's completed, activate the environment
conda activate AMR++_env

# You now have access to all the AMR++ software dependencies (locally)
samtools --help

# Run command to perform the demonstration pipeline using the conda profile.
nextflow run main_AMR++.nf


```
Now, you can check out the results in the newly created "test_results" directory.

# Using AMR++ to analyze your data

AMR++ is customizable to suit your computing needs and analyze your data. Primarily, the ```-profile``` paramater allows you to choose between running AMR++ using a singularity container, docker container, anaconda packages, or a local installation of your software. 
All parameters used to control how AMR++ analyzes your data can also be changed as needed in a variety of ways. For full information, review this [configuration document.](docs/configuration.md)


Below is a brief example, the default parameters were run using this command (with the conda environment, AMR++_env, already activated):

```nextflow run main_AMR++.nf```

To change the reads that were analyzed, you should specify the ```--reads`` parameters. Here, we can use regular expressions to point to your samples in a different directory.
```bash
nextflow run main_AMR++.nf --reads "path/to/your/reads/*_R{1,2}.fastq.gz" 
```


## Choosing the right pipeline

AMR++ analyzes data by combining workflows that takes a set of sequencing reads through various bioinformatic software. We recommend our standard AMR++ pipeline as a comprehensive way to start from raw sequencing reads, QC assessment, host DNA removal, and resistome analysis with MEGARes. However, users might only want to replicate portions of the pipeline and have more control over their computing needs. Using the ```--pipeline``` parameter, users can now change how AMR++ runs.

Check out [this document](docs/choosing_pipeline.md) for more details and guidance on picking the right ```--pipeline``` parameter.

## Running analyses in steps

Realistically, running the entire pipeline can be challenging due to storage limitations. Instead, we recommend running the pipeline in steps, which allows for erasing the "work" directory in between analytic steps. Remember, the work directory is only needed in case the pipeline run fails and you want to use ```-resume``` to pick up where you left off. 

Here are some tutorials to run each analysis step by step:

- [Paired-end analysis step-by-step](docs/Step_by_step_tutorial.md)
- [Single-end analysis step-by-step](docs/SingleEnd_read_tutorial.md)
- [Merged-read analysis step-by-step](docs/Merged_read_tutorial.md)

> **Tip:** Always run the demo first (`--pipeline demo`) before your first analysis. See [Getting Started](docs/GettingStarted.md) for an explanation of the work directory, SLURM setup, and conda vs local profiles.


# Optional flags

## SNP verification

AMR++ now works in conjuction with a [custom SNP verification software](https://github.com/Isabella136/AmrPlusPlus_SNP) to evaluate alignments to gene accessions requiring SNP confirmation to confer resistance. To include this workflow, include the ```--snp Y``` flag in your command like this:

```bash
nextflow run main_AMR++.nf -profile conda --snp Y
```
This will create with the standard count table (AMR_analytic_matrix.csv) in addition to a count matrix with SNP confirmed counts (SNPconfirmed_AMR_analytic_matrix.csv).

## Deduplicated counts

Another option is to include results for deduplicated counts by using the ```--deduped Y``` flag in your command.

```bash
nextflow run main_AMR++.nf -profile conda --snp Y --deduped Y
```

With this flag, AMR++ will extract the deduplicated alignments to MEGARes also output a count matrix with deduplicated counts. Since also we included the ```--snp Y``` flag, we will end up with 4 total output count matrices.

## Rarefaction analyzer

The final optional analyis is to perform rarefaction on resistome counts to evaluate sequencing depth at all resistome annotation levels (i.e. Type, Class, Mechanism, Group, Gene). You can run this analysis by adding ```--rarefaction Y``` to your command or modifying the params.config file. 

```bash
nextflow run main_AMR++.nf -profile conda --snp Y --deduped Y --rarefaction Y
```

With rarefaction analysis, we'll create various figures to summarize sequencing depth and output figures in the "ResistomeAnalysis/Rarefaction/Figures" directory. 
