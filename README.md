Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

## DRAFT ## DRAFT

## AMR++ v3 updates

Updates
- SNP confirmation software
- pipeline structure
 - switch to nextflow DSL2 and collection of "pipeline" options
 - snakemake approach for optimizing storage requirements
- Includes updates to MEGARes v3



Brief tutorial for nextflow pipeline test run


Brief tutorial for snakemake pipeline test run


```


# Microbial Ecology Group (MEG)
(https://megares.meglab.org/)

Our international multidisciplinary group of scientists and educators is addressing the issues of antimicrobial resistance (AMR) and microbial ecology in agriculture through research, outreach, and education. By characterizing risks related to AMR and microbial ecology, our center will identify agricultural production practices that are harmful and can be avoided, while also identifying and promoting production practices and interventions that are beneficial or do no harm to the ecosystem or public health. This will allow society to realize “sustainable intensification” of agriculture.

# MEGARes and the AMR++ bioinformatic pipeline
(http://megares.meglab.org/amrplusplus/latest/html/v2/)

The MEGARes database contains sequence data for approximately 9,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AmrPlusPlus can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth.

With AMR++, you will obtain alignment count files for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

More Information
----------------

- [Installation](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/installation.md)
- [Usage](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/usage.md)
- [Configuration](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/configuration.md)
- [Accessing AMR++](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/accessing_AMR++.md)
- [Output](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/output.md)
- [Dependencies](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/dependencies.md)
- [Software Requirements](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/requirements.md)
- [FAQs](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/FAQs.md)
- [Details on AMR++ updates](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/update_details.md)
- [Contact](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/contact.md)
