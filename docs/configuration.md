Configuration
-------------

The pipeline source code comes with a configuration file that can be used to set environment variables or default command-line options. Setting these variables before hand may be useful in situations when you do not want to specify a long list of options from the command line. This configuration file can be found in the root source code directory and is called **nextflow.config**. You can modify this file, save the changes, and run the pipeline directly.


Customize Environment Variables
-------------------------------

By default, the pipeline uses the default minikraken database (~4GB) to classify and assign taxonomic labels to your sequences. As Kraken loads this database into memory, this mini database is particularly useful for people who do not have access to large memory servers. If you would like to use a custom database or the standard Kraken database (~160GB), you will need to build it yourself and modify the **MINIKRAKENDB** environment variable to point to its location on your machine.

```bash
env {
    /* Location of minikraken database */
    MINIKRAKENDB = "/opt/minikraken"

    /* Location of trimmomatic jar file */
    TRIMMOMATIC = "/opt/trimmomatic/Trimmomatic-0.36"
}
```

Customize Command-line Options
------------------------------

The params section allows you to set the different commmand-line options that can be used within the pipeline. Here, you can specify input/output options, trimming options, and algorithm options.

If you intend to run multiple samples in parallel, you must specify a glob pattern for your sequence data as shown for the **reads** parameter. For more information on globs, please see this related [article](https://en.wikipedia.org/wiki/Glob_(programming)).

```bash
params {
    /* Location of forward and reverse read pairs */
    reads = "data/raw/*_{1,2}.fastq"

    /* Location of host genome index files */
    host_index = ""

    /* Location of host genome */
    host = "data/host/chr21.fasta"

    /* Location of amr index files */
    amr_index = ""

    /* Location of antimicrobial resistance (AMR) database */
    amr = "data/amr/megares_database_v1.01.fasta"

    /* Location of amr annotation file */
    annotation = "data/amr/megares_annotations_v1.01.csv"

    /* Location of adapter sequences */
    adapters = "data/adapters/nextera.fa"

    /* Location of tab delimited adapter sequences */
    fqc_adapters = "data/adapters/nextera.tab"

    /* Output directory */
    output = "./test"

    /* Number of threads */
    threads = 16

    /* Trimmomatic trimming parameters */
    leading = 3
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
