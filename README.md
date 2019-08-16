Overview
--------

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AmrPlusPlus can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth. AmrPlusPlus can perform this analysis as well.

With microbiome analysis of shotgun metagenomic reads, you may want to focus on the presence of a particular species. To improve the accuracy of species level classification, the "main_kraken.nf" can be run on the non-host reads in the BAMtoFASTQ/ output from "main_amr_plus_plus.nf". This script will run kraken with and without the kraken filter, pull out the reads classified as the species of interest, and use blastn for secondary classification. The targeted species can be changed using the "species" flag, with the default being Salmonella enterica. 

With AmrPlusPlus, you will obtain count files for each sample that can be combined into a count matrix and analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

## Example commands:
nextflow main.nf --reads "PATH/to/files/*_R{1,2}_001.fastq.gz " --kraken_db $krakendir --output $outputDir --threads 1 -w $workingDIR --host $HOSTFASTA -profile local -species "Salmonella enterica"

/s/angus/index/common/tools/nextflow run main_confirm_SNP_counts.nf --sams "/s/angus/index/common/tools/AlignToAMR/*_concatenated.amr.alignment.sam" --output /s/angus/index/common/tools/SNP_confirmation --threads 1 -w /s/angus/index/common/tools/work --host /s/angus/index/databases/bwa_indexes/mod_bos_taurus/mod_bos_taurus.fna -profile local -resume

/s/angus/index/common/tools/nextflow run main_AMR_amrplusplus_snp.nf --reads "/media/AngusWorkspace/MEGmobile_test/test_run/*_{1,2}.fastq" --output /media/AngusWorkspace/Megmobile_test --threads 10 -w /media/AngusWorkspace/work_mobile --host /s/angus/index/databases/bwa_indexes/mod_bos_taurus/mod_bos_taurus.fna --amr /s/angus/index/databases/MEGMobile/MEGmobile_database_10Feb2019.fasta -profile local_angus -resume


More Information
----------------

- [Software Requirements](https://github.com/cdeanj/amrplusplus/blob/master/docs/requirements.md)
- [Installation](https://github.com/cdeanj/amrplusplus/blob/master/docs/installation.md)
- [Usage](https://github.com/cdeanj/amrplusplus/blob/master/docs/usage.md)
- [Configuration](https://github.com/cdeanj/amrplusplus/blob/master/docs/configuration.md)
- [Output](https://github.com/cdeanj/amrplusplus/blob/master/docs/output.md)
- [Dependencies](https://github.com/cdeanj/amrplusplus/blob/master/docs/dependencies.md)
- [Contact](https://github.com/cdeanj/amrplusplus/blob/master/docs/contact.md)
