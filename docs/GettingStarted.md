
Getting started with AMR++
-----------------

To get started, view the [Installation document](https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/installation.md) to determine the best way to install AMR++ to your computing cluster.

Next, we will run a small sample dataset that comes with the pipeline source code. As such, we will not be specifying any input paths as they have already been included. In this example, we'll assume that either all tools are installed and in your $PATH (or you installed the conda environment and activated the environment).


```bash
# If you followed the instructions on the installation document, you must now navigate to the AMR++ directory
cd AMRplusplus

# Run test of AMR++ by specifying the "-profile local".
nextflow run main_AMR++.nf -profile local

# Explore the "test_results/" directory to view pipeline outputs
ls test_results/

```


We can now change the "--pipeline" parameter to perform a different set of analyses. View the [configuration doc](#https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/configuration.md) for more details on pipeline options and modifying further parameters. 

```bash
# Run another test, but now we add "--pipeline eval_qc" to only calculate QC stats.
# We'll also change the "--output" flag to store the results in a new directory.
nextflow run main_AMR++.nf -profile local --pipeline eval_qc --output test_QC_stats

# View the results
ls test_QC_stats/

# Look at the multiqc_report.html and modify the params.config file to change trimming parameters.

# Now, you can run just the QC trimming step with trimmomatic, by using the "--pipeline trim_qc" flag. Change output and work directory.
nextflow run main_AMR++.nf -profile local --pipeline trim_qc --output test_QC_trimming -w work_trim

# If the pipeline completes without issue, remember to erase your "work" directories to save storage space.

```


Alternatively, you can run the entire pipeline as shown below.

```bash

# You can use the "--pipeline standard_AMR" to run the standard AMR++ pipeline.
nextflow run main_AMR++.nf -profile local --pipeline standard_AMR --output test_AMR++_output -w work_AMR++

```
