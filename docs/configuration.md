# Contents

* [Configuration](#configuration)
* [Customize environmental variables using profiles](#customize-environment-variables-using-profiles)
* [Customize parameters using the commandline](#customize-amr-pipeline-parameters)
  * [Modifying the params.config file](#modifying-the-paramsconfig-file) 
  * [Modifying parameters using the command-line](#modifying-parameters-using-the-command-line)
    * [Analyzing your samples](#analyzing-your-samples)
    * [Running with Kraken](#running-with-kraken)
    * [Including SNP confirmation](#running-with-snp-confirmation)
    * [Including deduplicated count results](#running-with-deduplicated-counts)
* [Selecting the right pipeline](#selecting-the-right-pipeline)

## Configuration
-------------

The pipeline source code comes with two configuration files that can be used to set environment variables and default command-line options. These configuration files can be found in the root source code directory and are called **nextflow.config** and **params.config**.

The **nextflow.config** file mainly contains parameters regarding how AMR++ will run on your computing cluster using the ```--profile``` parameter. 

The **params.config** contains parameters that control which files are being analyzed and parameters for the software in the pipeline. Setting the variables in the **params.config** before hand may be useful in situations when you do not want to specify a long list of options from the command line or want to have a seperate file for each project. You can modify these files, save the changes, and run the pipeline directly. More details below.


## Customize Environment Variables using profiles
----------------------------------------------

The **nextflow.config** contains a section that allows the use of environment "profiles" when running AmrPlusPlus. Further information for each profile can be found within the /config directory. In brief, profiles allow control over how the pipeline is run on different computing clusters. We recommend the "singularity" profile which employs singularity containers which contain all the required bioinformatic tools.

We make the following profiles available to suit your computing needs; "local", "local_slurm", "conda","conda_slurm", "singularity", "apptainer", "singularity_slurm", and "docker". You specify which profile to use with the ```-profile`` flag.


```bash
profiles {
  local {
    includeConfig "config/local.config"
  }
  local_slurm {
    includeConfig "config/local_slurm.config"
    process.executor = 'slurm'
  }
  conda {
    includeConfig "config/conda.config"
    conda.enabled = true
    conda.cacheDir = "$baseDir/envs/"
    conda.useMamba = true
    conda.createTimeout = '30 min'
  }
  docker {
    includeConfig "config/local.config"
    docker.enabled = true
    process.container = 'enriquedoster/amrplusplus:latest'
  }
  singularity {
    includeConfig "config/singularity.config"
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.cacheDir = "$baseDir/envs/"
  }
   apptainer {
    includeConfig "config/apptainer.config"
    apptainer.enabled = true
    apptainer.autoMounts = true
    apptainer.cacheDir = "$baseDir/envs/"
  }
  conda_slurm {
    includeConfig "config/conda_slurm.config"
    process.executor = 'slurm'
    conda.cacheDir = "$baseDir/envs/"
    conda.enabled = true
    conda.useMamba = true
    conda.createTimeout = '30 min'
  }
   singularity_slurm {
    includeConfig "config/singularity_slurm.config"
    process.executor = 'slurm'
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.cacheDir = "$baseDir/envs/"
  }
}
```

## Customize AMR++ pipeline parameters
------------------------------

The params section allows you to set the different commmand-line options that can be used within the pipeline. Here, you can specify input/output options, trimming options, and algorithm options.

### Modifying the params.config file
Below is a list of all of the parameters that AMR++ uses by default. They can be found in the ```params.config``` file in the main directory. These parameters can be modified by changing this file or specifying any of these parameters on the command line using a double dash, like this: ```--reads "path/to/your/reads/*_R{1,2}.fastq.gz"```. Otherwise, change the parameters in the ```params.config``` file prior to running the AMR++ pipeline.

These are all of the parameters used by AMR++:
```bash
params {
    /* Location of forward and reverse read pairs */
    reads = "${baseDir}/data/raw/*_R{1,2}.fastq.gz"

    /* Location of reference/host genome */
    reference = "${baseDir}/data/host/chr21.fasta.gz"

    /* Output directory */
    output = "test_results"
    
    /* Kraken database location, default is "null" */   
    kraken_db = null

    /* Location of amr index files */
    amr_index = ""

    /* Location of antimicrobial resistance (MEGARes) database */
    amr = "${baseDir}/data/amr/megares_database_v3.00.fasta"

    /* Location of amr annotation file */
    annotation = "${baseDir}/data/amr/megares_annotations_v3.00.csv"

    /* Location of SNP confirmation script */
    snp_confirmation = "${baseDir}/bin/snp_confirmation.py"

    /* Number of threads */
    threads = 4

    /* Trimmomatic trimming parameters */
    adapters = "${baseDir}/data/adapters/nextera.fa"

    leading = 3
    trailing = 3
    slidingwindow = "4:15"
    minlen = 36

    /* Resistome threshold */
    threshold = 10

    /* Starting rarefaction level */
    min = 5

    /* Ending rarefaction level */
    max = 100

    /* Number of levels to skip */
    skip = 5

    /* Number of iterations to sample at */
    samples = 1

    /* multiQC */
    multiqc = "$baseDir/data/multiqc"

    /* Display help message */
    help = false
}
```
### Modifying parameters using the command-line

#### Analyzing your samples
------
If you intend to run multiple samples in parallel, you must specify a glob pattern for your sequence data as shown for the **reads** parameter. For more information on globs, please see this related [article](https://en.wikipedia.org/wiki/Glob_(programming)).

For example, the default parameters can be used to run the pipeline with this command:

```bash
nextflow run main_AMR++.nf -profile singularity
```

This will run the default samples through the pipeline and this can be seen below, under the ```--reads``` parameter. To change the reads that were analyzed, you should specify the ```--reads`` parameter on the command line. Here, we can use regular expressions to point to your samples in a different directory.

```bash
nextflow run main_AMR++.nf -profile singularity  --reads "path/to/your/reads/*_R{1,2}.fastq.gz" 
```

#### Running with Kraken
-----
By default, the pipeline uses the default minikraken database (~4GB) to classify and assign taxonomic labels to your sequences. As Kraken loads this database into memory, this mini database is particularly useful for people who do not have access to large memory servers. We provide a script to easily download the minikraken database.

```bash
 sh download_minikraken.sh
```

If you would like to use a custom database or the standard Kraken database (~160GB), you will need to build it yourself and modify the **kraken_db** environment variable in the ```params.config ``` file to point to its location on your machine. 

#### Running with SNP confirmation
-----
To include SNP confirmation as part of the AMR++ analysis, you have to include the ```--snp Y``` flag. Like this:

```bash
nextflow run main_AMR++.nf -profile singularity  --reads "path/to/your/reads/*_R{1,2}.fastq.gz" --snp Y
```

#### Running with deduplicated counts
-----
Additionally, you can also output deduplicated counts by cinluding the flag, ```--deduped Y```. Like this:

```bash
nextflow run main_AMR++.nf -profile singularity  --reads "path/to/your/reads/*_R{1,2}.fastq.gz" --snp Y --deduped Y
```


## Selecting the right pipeline

AMR++ now includes the option to run different components of the pipeline at a time by specifying the ```--pipeline``` flag.

Main pipeline options
  * Standard AMR pipeline ( QC trimming > Host DNA removal > Resistome alignment > Resistome results)
    ```bash
    --pipeline standard_AMR
    ```
  * Fast AMR pipeline (QC trimming > Resistome alignment > Resistome results)
    ```bash
    --pipeline fast_AMR
    ```
  * AMR pipeline with kraken ( QC trimming > Host DNA removal > Resistome alignment > Resistome results) & (Non-host reads > Microbiome analysis)
    ```bash
    --pipeline standard_AMR_wKraken
    ```
  * 16S Microbiome analysis with qiime2 (DADA2 QC > Classification with SILVA)
    ```bash
    --pipeline qiime2
    ```
    Pipeline components
  * Evaluate QC with multiQC
    ```bash
    --pipeline eval_qc
    ```
  * QC trimming with trimmomatic
    ```bash
    --pipeline trim_qc
    ```
  * Align reads to host DNA and remove contaminants
    ```bash
    --pipeline rm_host
    ```
  * Only perform AMR++ resistome analysis
    ```bash
    --pipeline resistome
    ```
  * Only perform microbiome analysis with Kraken
    ```bash
    --pipeline kraken
    ```