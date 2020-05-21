Details on AMR++ updates
------------

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
