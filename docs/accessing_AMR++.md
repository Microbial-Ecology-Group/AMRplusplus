Accessing AMR++
------------

This section will help you get access to all the bioinformatic tools required for metagenomic analysis with AMR++.

Amazon Web Services
-----

In order to facilitate evaluation of the MEGARes 2.0 database and the functionality of AMR++ 2.0 pipeline, we have provided free access to an Amazon Machine Image (AMI) with example files for analysis. AMR++ 2.0 is pre-installed and fully integrated with all necessary bioinformatic tools and dependencies within an AMI named "Microbial_Ecology_Group_AMR_AMI", allowing users to easily employ the AMR++ v2.0 pipeline within the Amazon Web Services (AWS) ecosystem. Please follow the instructions on amazon web services for details on creating your own EC2 instance (https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html). With this approach, users pay for the cost of a suitable AWS EC2 instance without the challenge of accessing large computing clusters and individually installing each piece of software necessary to run the pipeline (including all dependencies). Integration within AWS also allows users to scale the computing resources to fit the needs of any project size.

Singularity container
-----------------

Singularity containers allow the packaging of multiple bioinformatic tools. While singularity is a popular tool and likely to be supported by many computing clusters, please contact your system administrator for help with installing singularity. Installation on a local computer is also an option and can be performed by following these instructions: https://sylabs.io/guides/3.0/user-guide/installation.html

We provide AMR++ with a singularity container that is automatically accessed when running the AMR++ pipeline by using the flag, "-profile singularity". Additionally, the singularity container is supported on singularity-hub.org and can be used locally for custom analysis (https://singularity-hub.org/collections/3418). P

```bash
# Choose your preference to pull the container from Singularity Hub (once)
$ singularity pull shub://meglab-metagenomics/amrplusplus_v2

# Then interact with it (enter "exit" to leave the singularity container):
$ singularity shell amrplusplus_v2.sif

```


