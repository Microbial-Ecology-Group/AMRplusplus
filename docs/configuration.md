Configuration
-------------

The pipeline source code comes with a configuration file that can be used to set environment variables or default command-line options. Setting these variables before hand may be useful in situations when you do not want to specify a long list of options from the command line. This configuration file can be found in the root source code directory and is called **nextflow.config**. You can modify this file, save the changes, and run the pipeline directly.


Customize Environment Variables using profiles
----------------------------------------------

The **nextflow.config** contains a section that allows the use of environment "profiles" when running AmrPlusPlus. Further information for each profile can be found within the /config directory. In brief, profiles allow control over how the pipeline is run on different computing clusters. We recommend the "singularity" profile which employs a Singularity container  with all the required bioinformatic tools.


```bash
profiles {
  local {
    includeConfig "config/local.config"
  }
  local_angus {
    includeConfig "config/local_angus.config"
  }
  local_MSI {
    includeConfig "config/local_MSI.config"
  }
  slurm {
    process.executor = 'slurm'
    includeConfig "config/slurm.config"
    process.container = 'shub://meglab-metagenomics/amrplusplus_v2'
  }
  singularity {
    includeConfig "config/singularity.config"
    process.container = 'shub://meglab-metagenomics/amrplusplus_v2'
  }
}
```

Customize Command-line Options
------------------------------

The params section allows you to set the different commmand-line options that can be used within the pipeline. Here, you can specify input/output options, trimming options, and algorithm options.

If you intend to run multiple samples in parallel, you must specify a glob pattern for your sequence data as shown for the **reads** parameter. For more information on globs, please see this related [article](https://en.wikipedia.org/wiki/Glob_(programming)).


By default, the pipeline uses the default minikraken database (~4GB) to classify and assign taxonomic labels to your sequences. As Kraken loads this database into memory, this mini database is particularly useful for people who do not have access to large memory servers. We provide a script to easily download the minikraken database.

> sh download_minikraken.sh

If you would like to use a custom database or the standard Kraken database (~160GB), you will need to build it yourself and modify the **kraken_db** environment variable in the nextflow.config file to point to its location on your machine. 


```bash
params {
    /* Location of forward and reverse read pairs */
    reads = "data/raw/*_{1,2}.fastq.gz"

    /* Location of adapter sequences */
    adapters = "data/adapters/nextera.fa"

    /* Location of tab delimited adapter sequences */
    fqc_adapters = "data/adapters/nextera.tab"

    /* Location of host genome index files */
    host_index = ""

    /* Location of host genome */
    host = "data/host/chr21.fasta.gz"
    
    /* Kraken database location, default is "none" */   
    kraken_db = "minikraken2_v2_8GB_201904_UPDATE"

    /* Location of amr index files */
    amr_index = ""

    /* Location of antimicrobial resistance (MEGARes) database */
    amr = "data/amr/megares_database_v1.02.fasta"

    /* Location of amr annotation file */
    annotation = "data/amr/megares_annotations_v1.02.csv"

    /* Location of SNP metadata */
    snp_annotation = "data/amr/snp_location_metadata.csv"

    /* Location of SNP confirmation script */
    snp_confirmation = "bin/snp_confirmation.py"

    /* Output directory */
    output = "test_results"

    /* Number of threads */
    threads = 10
    smem_threads = 12

    /* Trimmomatic trimming parameters */
    leading = 10
    trailing = 3
    slidingwindow = "4:15"
    minlen = 36

    /* Resistome threshold */
    threshold = 80

    /* Starting rarefaction level */
    min = 5

    /* Ending rarefaction level */
    max = 100

    /* Number of levels to skip */
    skip = 5

    /* Number of iterations to sample at */
    samples = 1

    /* Display help message */
    help = false
}
```
