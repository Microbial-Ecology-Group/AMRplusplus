Details on AMR++ updates
------------

## 2026-01-19: AMR++ v3 update
- Added functionality for single end reads and a pipeline to merge reads with FLASH and analyze both merged fragments and unmerged reads.
- Resistome analysis: Changed default gene fraction threshold to "0" and changed extraction of alignments from the BAM file to only include primary alignments (one read = one alignment)
  - Based on evaluation of various datasets, we identified that secondary and supplemental alignments make up less than 2% of all alignments and only a small subset of reads contain more than a primary alignment.
  - Therefore we decided to simplify counting of alignments to only include primary alignments
- Updated slurm configuration profiles to include custom SBATCH parameters corresponding to each module in AMR++. This allows for submission of an SBATCH script with minimal resources to run AMR++ and each job will automatically be submitted with the required resources.
- Improved Kraken2 database handling and output processing
  - Fixed kraken2 long-to-wide format conversion
  - Added support for paired-end flag in kraken workflow
  - Added merged kraken workflow for analyzing merged reads
- Enhanced SNP verification workflow
  - Updated SNP results parsing for better handling of detailed output
  - Improved SNP confirmation code integration
  - Store detailed SNP verification output
- Improved documentation
  - Added comprehensive step-by-step tutorials
  - Created Single-End read tutorial
  - Updated installation documentation for Apptainer
  - Corrected parameter documentation (--host instead of --reference)
- Configuration improvements
  - Updated MultiQC integration and conda environment
  - Enhanced resource allocation defaults
  - Added QIIME workflow integration
  - Improved slurm profile labels and configuration
- Bug fixes and optimizations
  - Fixed host removal command and samtools flags
  - Improved sample ID naming and handling
  - Fixed kraken output saving and database selection
  - Cleaned up output directory structure
  - Fixed deduplication processes
  - Updated trimming parameters (added crop_len option)
  - Various fixes to process dependencies and naming conventions
- Maintenance updates
  - Updated to Python 3.9 in Docker container
  - Cleaned up container documentation
  - Updated help messages and usage documentation 


## 2022-09-06 : AMR++ v3 update
- Change in repository from [AMR++ v2](https://github.com/meglab-metagenomics/amrplusplus_v2) to this repository under the [microbial ecology group github page](https://github.com/Microbial-Ecology-Group/AMRplusplus). This repository will include all further updates to AMR++. 
- Addition of [SNP confirmation software](https://github.com/Isabella136/AmrPlusPlus_SNP)
- pipeline structure
   - switch to nextflow DSL2 and collection of "pipeline" options
   - snakemake approach for optimizing storage requirements
- Includes updates to MEGARes v3


## 2020-05-21 : AMR++ v2.0.2 update
Fixed a mistake with the config/singularity.config file to correctly call the singularity container anytime that RGI is run.

## 2020-05-21 : AMR++ v2.0.1 update
We identified issues with running RGI thanks to github users, AroArz and DiegoBrambilla. As of this update, RGI developers are focused on contributing to the COVID-19 response, so we plan to reconvene with them when their schedule opens up. In the meantime, we are releasing updates to continue AMR++ functionality.
We found that the errors were associated with RGI bugs that were previously reported:
* In [RGI issue #93](https://github.com/arpcard/rgi/issues/93), the github user, mahesh-panchal, reported that you need to run the "rgi main" command twice on the same dataset for it to successfully complete the analysis. The RGI component of AMR++ has been updated to work for now, but we plan for further changes to clean up the code.
* In [RGI issue #60](https://github.com/arpcard/rgi/issues/60), caspargross reported issues with containerizing RGI due to requirements for a "writable file system". As a temporary fix, we updated AMR++ code so that the user has to download the CARD database locally and use an additional flag to specify the location of the local database.
* Errors in running RGI will now be "ignored" so that the pipeline continues running but still provides any temporary files created with RGI. This should allow you to troubleshoot on your own or run any additional analysis using the reads aligning to gene accessions that "RequireSNPConfirmation"

Other updates:
* minor fixes to singularity/slurm configuration files
* updated the "resistome" and "rarefaction" code. It is now included in the "bin" directory.
* updated the script for downloading the latest version of mini-kraken
* created new singularity container just for the RGI software
* output zipped nonhost files directly
