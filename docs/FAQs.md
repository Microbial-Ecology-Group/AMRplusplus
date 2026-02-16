# Troubleshooting and Frequently Asked Questions (FAQs)

This document covers common issues you may encounter when running AMR++ and provides guidance for diagnosing and resolving errors.

## Table of Contents

- [Basic Troubleshooting](#basic-troubleshooting)
  - [Reading the Nextflow Log](#reading-the-nextflow-log)
  - [Inspecting the Work Directory](#inspecting-the-work-directory)
- [Common Errors](#common-errors)
  - [Tool Not Found / Command Not Found](#tool-not-found--command-not-found)
  - [Dependency Download Fails (No Internet on Compute Node)](#dependency-download-fails-no-internet-on-compute-node)
  - [Out of Memory Errors](#out-of-memory-errors)
  - [File Permission Errors](#file-permission-errors)
- [General FAQs](#general-faqs)
  - [Which Profile Should I Use?](#which-profile-should-i-use)

---

## Basic Troubleshooting

When a pipeline run fails, Nextflow provides several ways to figure out what went wrong. The two most important things to check are the **Nextflow log file** and the **work directory** for the failed task.

### Reading the Nextflow Log

Every Nextflow run generates a hidden log file called `.nextflow.log` in the directory where you launched the pipeline. This file contains detailed information about every step of the run, including the exact error messages and the work directory paths for each task.

To find the error, search for `ERROR` or the name of the failed process:

```bash
grep -i "ERROR" .nextflow.log
```

You can also look for the work directory of the failed task:

```bash
grep "work/" .nextflow.log | tail -20
```

The log will show lines like:

```
[ab/cd1234] process > runresistome (SampleA) [100%] 1 of 1, failed: 1
```

The hash in brackets (`ab/cd1234`) corresponds to a subdirectory under `work/` where that specific task ran. This is where you go to investigate further.

### Inspecting the Work Directory

Each task that Nextflow executes runs inside its own directory under `work/`. For the example above, the full path would be something like:

```
work/ab/cd1234e56789...
```

You can find the full path in `.nextflow.log`, or use `nextflow log` to list recent runs and their task directories. Inside each work directory, you will find:

| File | Description |
|------|-------------|
| `.command.sh` | The exact shell script that Nextflow ran for this task. This is useful for understanding what command was executed and for manually re-running it to reproduce the error. |
| `.command.run` | The wrapper script Nextflow uses to execute the task (includes environment setup). |
| `.command.out` | Standard output (stdout) from the task. |
| `.command.err` | Standard error (stderr) from the task — **check this first** for error messages. |
| `.command.log` | Combined log output (used by some executors). |
| `.exitcode` | The exit code of the task. `0` means success; anything else indicates a failure. |

A typical debugging workflow looks like this:

```bash
# 1. Find the work directory for the failed task
grep "failed" .nextflow.log | tail -5

# 2. Navigate to it
cd work/ab/cd1234e56789

# 3. Check the error output
cat .command.err

# 4. Look at the exact command that was run
cat .command.sh

# 5. Check the exit code
cat .exitcode
```

If you want to re-run the failed command manually to test a fix, you can execute `.command.sh` directly from within the work directory. Just make sure the necessary tools are available in your environment.

---

## Common Errors

### Tool Not Found / Command Not Found

**Symptom:** The error log shows something like:

```
.command.sh: line 3: bwa: command not found
```

or:

```
FATAL: cannot execute binary file
```

**Cause:** The bioinformatic tool required by the process is not installed or not available in your `$PATH`.

**Solution:** We recommend setting up the conda environment before submitting your pipeline job. This ensures all tools are installed and available:

```bash
# 1. Create the environment first (only needed once)
mamba env create -f envs/AMR++_env.yaml

# 2. Activate it
conda activate AMRplusplus

# 3. Verify key tools are available
bwa 2>&1 | head -3
samtools --version | head -1

# 4. Now submit your pipeline job (the environment carries over)
sbatch run_AMR++_slurm.sbatch
```

If you are using the `conda` profile, make sure `mamba` is available and that the conda cache directory (`envs/`) is writable. On some clusters, the first run may take a while as the environment is built. Subsequent runs will reuse the cached environment.

For Singularity or Apptainer profiles, the container image includes all necessary tools — you should not need to install anything separately.

### Dependency Download Fails (No Internet on Compute Node)

**Symptom:** The `build_dependencies` process fails with errors like:

```
fatal: unable to access 'https://github.com/Isabella136/AmrPlusPlus_SNP.git/': Could not resolve host: github.com
```

**Cause:** The `build_dependencies` process needs to clone a Git repository (the SNP verification tool), but compute nodes on many HPC clusters do not have internet access. This commonly happens when you submit an `sbatch` job and the job lands on a compute node that is isolated from the network.

**Solutions:**

1. **Run the demo first from a login node** (recommended). Login nodes typically have internet access. Running the demo pipeline once will download and cache the dependencies locally so that subsequent runs on compute nodes can find them:

   ```bash
   # On the login node (with internet access):
   nextflow run main_AMR++.nf -profile conda --pipeline demo
   ```

   After this completes, the SNP verification files will be in `bin/AmrPlusPlus_SNP/` and future runs will skip the download.

2. **Load a web proxy module.** Some clusters provide a module that enables internet access on compute nodes. Check with your HPC administrators — the module name varies by system, but it is often something like:

   ```bash
   module load WebProxy
   ```

   Add this line to your sbatch script *before* the `nextflow run` command.

3. **Clone the repository manually.** If neither option above works, you can download the dependency yourself from a node with internet access:

   ```bash
   cd bin/
   git clone https://github.com/Isabella136/AmrPlusPlus_SNP.git
   chmod -R 777 AmrPlusPlus_SNP/
   ```

### Out of Memory Errors

**Symptom:** A process fails with exit code 137 or 140, or the `.command.err` file shows:

```
slurmstepd: error: Detected 1 oom-kill event(s)
```

**Cause:** The task exceeded its allocated memory. This is most common with Kraken2 when using a large database, or with alignment steps on very large FASTQ files.

**Solution:** If you are using a `_slurm` profile, the pipeline is already configured to automatically retry failed tasks with 1.5× the memory. If the retry also fails, you can increase the memory allocation for that label in `config/local_slurm.config`. For example:

```groovy
withLabel: xlarge {
    memory = { 350.GB * (task.attempt > 1 ? 1.5 : 1) }
}
```

For Kraken2 specifically, you can also add `--memory-mapping` to `params.kraken_options` in `params.config`. This drastically reduces memory usage by reading the database from disk instead of loading it into RAM, at the cost of significantly longer runtimes.

### File Permission Errors

**Symptom:** Errors like `Permission denied` when writing output files or accessing tools installed by other users.

**Cause:** On shared servers with multiple users, certain directories may restrict write or execute permissions. Common cases include the pipeline output directory being owned by another user, or bioinformatic tools installed system-wide without execute permissions for your group.

**Solution:** Check the permissions of the relevant directory:

```bash
ls -la /path/to/directory
```

The permission string (e.g., `-rwxr-xr--`) shows read (r), write (w), and execute (x) permissions for the owner, group, and others. If you need write access to a shared directory, contact your system administrator or choose an output directory within your own home or scratch space:

```bash
nextflow run main_AMR++.nf --output /scratch/$USER/amr_results ...
```

For more information on Linux file permissions, see [this tutorial](https://www.guru99.com/file-permissions.html).

---

## General FAQs

### Which Profile Should I Use?

The correct profile depends on your computing environment:

| Environment | Recommended Profile |
|-------------|-------------------|
| Tools already installed and in `$PATH` | `local` |
| Tools installed, using SLURM | `local_slurm` |
| No tools installed, have conda/mamba | `conda` |
| Conda with SLURM | `conda_slurm` |
| Have Singularity installed | `singularity` or `singularity_slurm` |
| Have Apptainer installed | `apptainer` |
| Have Docker installed | `docker` |

If you have Singularity or Apptainer available on your system, these are often the easiest option since the container includes all dependencies and avoids installation issues entirely. If you are running many samples or using a large Kraken database, we recommend using a `_slurm` profile — see the [SLURM Resource Labels](Running_with_SLURM.md) documentation for details.

Also view this document for more [details on picking the right profile](choosing_pipeline.md).

