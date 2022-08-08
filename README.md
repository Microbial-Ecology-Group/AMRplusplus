Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)


## AMR++ v2 updates


### 2022-08-22 : AMR++ update coming soon
Hello AMR++ users, we would like to sincerely apologize for the delay in addresssing your concerns and updating AMR++. As a lot of you likely experienced, COVID was challenging and we were not able dedicate the resources to AMR++ that it deserves. We are happy to announce that we have assembled a team for another major update to AMR++ and the MEGARes database in the next few months!

A few notes:
 * We are aware of the issues with integrating RGI results with the AMR++ pipeline. Unfortunately, we are discontinuing our support of integrating AMR++ results with the RGI software.
 * We are attempting to remedy the issues that AMR++ users have reported, but we would also like to hear any other suggestions you might have. Please send any suggestions to enriquedoster@gmail.com with the subject line, "AMR++ update".
 * A few upcoming updates: easy control over the amount of intermediate files that are stored, option to re-arrange pipeline processes, better sample summary statistics provided, and improved functionality through nextflow profiles.


### 2020-03-21 : AMR++ v2.0.2 update.
We identified issues in running RGI with the full AMR++ pipeline thanks to github users, AroArz and DiegoBrambilla. We are releasing v2.0.1 to continue AMR++ functionality, but we are planning further updates for the next stable release. As of this update, RGI developers are focused on contributing to the COVID-19 response, so we plan to reconvene with them when their schedule opens up.
  * Please view the [CHANGELOG](https://github.com/meglab-metagenomics/amrplusplus_v2/blob/master/docs/CHANGELOG.md) for more details on changes included in AMR++ v2.0.1
  * To run the AMR++ pipeline with RGI, you'll have to download the CARD database locally and specify it's location using the "--card_db" flag like this:

```
# If you want to include RGI in your analysis, first download CARD with this command:
# We tested AMR++ v2.0.2 with the CARD database v3.0.8, but we recommend using the command below to get the latest CARD db
wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2

# In case the latest CARD database is causing issues, you can download the version we used for testing, v3.0.8:
wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/download/0/broadstreet-v3.0.8.tar.bz2 && tar xfvj card-data.tar.bz2


# If you run into an error regarding "Issued certificate has expired.", try this command:
wget --no-check-certificate -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2


# Run the AMR++ pipeline with the "--card_db" flag
nextflow run main_AmrPlusPlus_v2_withRGI.nf -profile singularity --card_db /path/to/card.json --reads '/path/to/reads/*R{1,2}_001.R1.fastq.gz' --output AMR++_results -w work_dir
```


# Microbial Ecology Group (MEG)
(https://megares.meglab.org/)

Our international multidisciplinary group of scientists and educators is addressing the issues of antimicrobial resistance (AMR) and microbial ecology in agriculture through research, outreach, and education. By characterizing risks related to AMR and microbial ecology, our center will identify agricultural production practices that are harmful and can be avoided, while also identifying and promoting production practices and interventions that are beneficial or do no harm to the ecosystem or public health. This will allow society to realize “sustainable intensification” of agriculture.

# MEGARes and the AMR++ bioinformatic pipeline
(http://megares.meglab.org/amrplusplus/latest/html/v2/)

The MEGARes database contains sequence data for approximately 8,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

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
