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
  }
  conda {
    includeConfig "config/conda.config"
    conda.enabled = true
    conda.cacheDir = "$baseDir/envs/"
    conda.useMamba = true
    conda.createTimeout = '30 min'
  }
  conda_slurm {
    includeConfig "config/conda_slurm.config"
    conda.cacheDir = "$baseDir/envs/"
    conda.enabled = true
    conda.useMamba = true
    conda.createTimeout = '30 min'
  }
  docker {
    includeConfig "config/docker.config"
    docker.enabled       = true
    docker.runOptions     = '-u $(id -u):$(id -g)'
    singularity.enabled  = false
    apptainer.enabled    = false
  }
  docker_slurm {
    includeConfig "config/docker_slurm.config"
    docker.enabled       = true
    docker.runOptions     = '-u $(id -u):$(id -g)'
    singularity.enabled  = false
    apptainer.enabled    = false
  }
  singularity {
    includeConfig "config/singularity.config"
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.cacheDir = "$baseDir/envs/"
  }
  singularity_slurm {
    includeConfig "config/singularity_slurm.config"
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
  apptainer_slurm {
    includeConfig "config/apptainer_slurm.config"
    apptainer.enabled = true
    apptainer.autoMounts = true
    apptainer.cacheDir = "$baseDir/envs/"
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
    /* Default pipeline */
    pipeline = "demo"

    /* Display help message */
    help = false

    /* Location of forward and reverse read pairs */
    reads = "${baseDir}/data/raw/*_R{1,2}.fastq.gz"

    /* for merged analysis */
    merged_reads = 'test_results/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz'

    /* Output directory */
    output = "test_results"

    /* Optional input for bam files for use with "--pipeline bam_resistome" */
    bam_files = null

    /* Default memory to run clumpify */
    clumpify_mem_gb = 8

    /* Location of reference/host genome */
    host = "${baseDir}/data/host/chr21.fasta.gz"

    /* Optionally, you can specify the location of the host index files created with bwa with the path and wildcard (*): */
    /* If you don't have the index files, replace this with "null" without quotes */
    host_index =  null
    
    /* Kraken database location, default is "null" */   
    kraken_db = null
    
    /* Kraken confidence score, 0.0 by default */
    kraken_confidence = 0.0

    kraken_options = ""

    /* Location of amr index files with wildcard */
    amr_index = "${baseDir}/data/amr/megares_database_v3.00.fasta*"

    /* Location of antimicrobial resistance (MEGARes) database */
    amr = "${baseDir}/data/amr/megares_database_v3.00.fasta"

    /* Location of amr annotation file */
    annotation = "${baseDir}/data/amr/megares_annotations_v3.00.csv"

    /* Samtools Resistome alignment flag options */
    samtools_flag = "" // For example: "-F 2304" to remove secondary and supplemental alignments

    /* Add SNP analysis */
    snp = "Y"

    /* Resistome threshold - associated with outdated resistomeanalyzer, now replaced with alignment_analyzer.py */
    threshold = 0

    /* Add rarefaction analysis */ 
    rarefaction = "N"

    /* Number of threads */
    threads = 4

    /* Trimmomatic trimming parameters */
    adapters = "${baseDir}/data/adapters/nextera.fa"
    leading = 3
    trailing = 3
    slidingwindow = "4:15"
    minlen = 36
    crop_len = 200

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

    /* Optional read deduplication (Clumpify) after QC trimming */
    read_dedup            = "N"         // Y inserts dedup between trim and host removal

    /* Add read deduplicaation analysis */
    deduped = "N"
    prefix = "AMR"

    /* Qiime2 - convenience subworkflow, not primary objective of AMR++ */
    /* Dada parameters */
    p_trim_left_f = 25
    p_trim_left_r = 26
    p_trunc_len_f = 225
    p_trunc_len_r = 220

    /* qiime2 bayes classifier */
    dada2_db = "$baseDir/data/qiime/gg-13-8-99-515-806-nb-classifier.qza"

    /* ── Alignment filtering (applied by alignment_analyzer.py) ───────────
     * Three independent filters, all on a 0 to 1 proportion scale:
     *
     *   min_gene_fraction    how much of the GENE must be covered
     *   min_query_coverage   how much of the READ must align
     *   min_identity         how well the ALIGNED PORTION must match
     *
     * match_qcov changes HOW min_query_coverage is calculated. The threshold
     * value means the same thing either way; only the metric changes:
     *   N  aligned_length / read_length            (mismatches count as covered)
     *   Y  (aligned_length - NM) / read_length     (genuine matches only)
     *
     * Y requires the NM tag. With min_query_coverage > 0 under match_qcov=Y,
     * or with any min_identity > 0, alignments lacking NM are excluded.
     */
    count_mode            = "fragment"  // DEFAULT. Mates resolved to one fragment.
                                        // Alternatives: read_end, alignment.
    group_aware           = "Y"         // DEFAULT ON. Same-Group mate disagreements -> 1 hit
                                        // (tie-break: higher match_qcov).
    edge_aware_qcov       = "Y"         // DEFAULT ON. Coverage-anchored qcov. 
    include_supplementary = "N"
    include_secondary = "N"
    cigar_aware_coverage  = "Y"
    per_read_alignment_stats = "N"      // Optional file with alignment data for each read. 

    /* Filters, all proportions 0 to 1 */
    min_gene_fraction     = 0
    min_query_coverage    = 0
    min_identity          = 0
    min_mapq              = 0

    /* Changes HOW min_query_coverage is calculated. Threshold value is unchanged.
     * N: aligned_length / read_length      Y: (aligned_length - NM) / read_length */
    match_qcov            = "N"     


    /* ── Coverage threshold evaluation ─────────────────────────────────────────
     * The evaluation always tests a two-dimensional grid: gene fraction against one
     * read-level filter. All values are comma-separated proportions from 0 to 1,
     * and "0" turns a filter off.
     *
     * AXIS 1 is always gene fraction.
     *
     * AXIS 2 is ONE of the three below. Set your choice to a list of values and
     * leave the other two at "0":
     *
     *   evaluation_query_coverage  how much of the READ aligned, measured as
     *                         aligned_length / read_length
     *   evaluation_match_qcov      how much of the READ aligned, measured as
     *                         (aligned_length - NM) / read_length, counting
     *                         only matching bases within the alignment
     *   evaluation_identity        how well the ALIGNED PORTION matched, measured
     *                         as (aligned_length - NM) / aligned_length
     *
     * evaluation_query_coverage and evaluation_match_qcov are two ways of measuring the
     * SAME property, so evaluating both is meaningless. Choose whichever
     * definition of query coverage you want to filter on. evaluation_identity
     * measures something different and is a genuine alternative axis.
     *
     * A single NON-ZERO value is a FIXED filter applied at every point in the
     * grid rather than a evaluated axis. So evaluation_identity = "0.9" holds identity
     * at 90% throughout while you evaluation gene fraction against a query-coverage
     * measure.
     */

    /* AXIS 1: always evaluated */
    evaluation_gene_fraction   = "0,0.1,0.25,0.5,0.8"                    // proportion, 0 to 1

    /* AXIS 2: evaluation ONE of these three; leave the others at "0",
     * or give one a single value to apply it as a fixed filter */
    evaluation_query_coverage  = "0,0.5,0.6,0.7,0.8,0.9,0.95"            // proportion, 0 to 1
    evaluation_match_qcov      = "0"                                     // proportion, 0 to 1
    evaluation_identity        = "0"                                     // proportion, 0 to 1

    /* Applied at every point in the grid regardless of what is evaluated */
    evaluation_edge_aware_qcov = "Y"   // coverage-anchored qcov: soft and hard clips
                                  // treated alike, and read bases overhanging a
                                  // short gene's edge are not penalized
    evaluation_exclude_snp     = "Y"   // Y drops RequiresSNPConfirmation genes entirely


    /* ── SNV calling (bam_snv) ────────────────────────────────────────────
     * The BAMs must be aligned to the SAME fasta given here. Use the MEGARes
     * REPRESENTATIVE database: both the NGLess filter and metaSNV require
     * uniquely-mapping reads, and the complete database's near-identical
     * accessions turn a large share of short reads into multi-mappers.
     */
    snv_reference          = "${baseDir}/data/amr/megares_database_v3.00.fasta"

    /* NGLess alignment filters. Applied before variant calling.
     * min_match_size is the aligned block size in bp; 45 is the floor used in
     * the benchmarking workflow, 100 is the stricter setting.
     * min_identity_pc is percent ANI; metaSNV guidance is 97 or above. */
    snv_min_match_size     = 100
    snv_min_identity_pc    = 97

    /* metaSNV options */
    snv_n_splits           = 1     // >1 writes one raw SNP file per split
    snv_db_ann             = ""    // optional gene annotation file
    snv_snpfile_prefix     = "called_SNPs"   // prefix shared by the raw call files

    /* Which alignment set was used; controls output filenames */
    snv_aln_wf             = "Standard"      // or "Deduped"

    /* Restrict the SNV matrix to genes also detected in the resistome matrix.
     * Requires a resistome run, so leave N when running bam_snv alone. */
    snv_filter_by_resistome = "N"
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

#### Running with deduplicated counts or reads
-----
You can run the AMR++ pipeline and have it deduplicate reads after read QC trimming and before host alignment and removal by including the flag ```read_dedup Y```.

```bash
nextflow run main_AMR++.nf -profile singularity  --reads "path/to/your/reads/*_R{1,2}.fastq.gz" --snp Y --read_dedup Y
```

Additionally, you can also output deduplicated alignments by including the flag, ```--deduped Y```. Like this:

```bash
nextflow run main_AMR++.nf -profile singularity  --reads "path/to/your/reads/*_R{1,2}.fastq.gz" --snp Y --deduped Y
```



## Selecting the right pipeline

AMR++ now includes the option to run different components of the pipeline at a time by specifying the ```--pipeline``` flag.

For more information about picking the right pipeline, [read this document](choosing_pipeline.md). Below is a brief summary of some options.

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