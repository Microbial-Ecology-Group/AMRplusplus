# Overview of steps in AMR++ with merged reads

# Confirm everything works

You'll need to load the typical AMR++ environment, plus add "FLASH" and "SeqKit". Or you can simply re-install the AMR++ conda environment using the recipe as shown in the installation tutorial. 

You'll still need to run the "demo" as this downloads the github with the SNP confirmation software. If you forget to do this, make sure your sbatch script includes the module that allows internet connection (e.g. "module load WebProxy").

## Download a new AMR++ repository and switch to "dev" branch

```
git clone https://github.com/Microbial-Ecology-Group/AMRplusplus

cd AMRplusplus/

git pull origin dev
```

# Full pipeline

You can run the full AMR++ pipeline (without kraken) by specifying the path to your `--reads` and using the `--pipeline se_AMR` flag. This could work well if you re-direct the results or the work directory `-w /path/to/shared_drive`.

## Step 1 : Run the "eval_qc" pipeline as normal with `--reads`

Example command:
```
nextflow run main_AMR++.nf --pipeline eval_qc --output SE_AMR++_analysis --reads "data/raw/*.fastq.gz" -profile local
``` 

## Step 2: Run "se_trim_qc" pipeline with `--reads`
Modify trimming parameters as needed in the `params.txt` file. 
Notice the `--reads` flag does not require a regular expression pattern to match paired reads, we can just use the asterisk "*" and file extension name to point to all samples. 

Example command:
```
nextflow run main_AMR++.nf --pipeline se_trim_qc --output SE_AMR++_analysis --reads "data/raw/*.fastq.gz" -profile local
``` 



## Step 3: Run "se_rm_host" with `--reads`

### Parameters to change

Parameters that have to change:
* `--pipeline` ==> `--pipeline rm_host`
* `--reads`  ==> `--reads "SE_AMR++_analysis/HostRemoval/NonHostFastq/*.fastq.gz"`
* `host` ==> `--host "/path/to/your/host/chr21.fasta.gz"` 
    * remember, you can change this in `params.config` file or add it to your nextflow command.
    * On grace, bovine: `/scratch/group/big_scratch/SHARED_resources/host_genome/GCF_002263795.3_ARS-UCD2.0_genomic.fna`

Example command:
```
nextflow run main_AMR++.nf --pipeline se_rm_host --output SE_AMR++_analysis --reads "SE_AMR++_analysis/HostRemoval/NonHostFastq/*.fastq.gz" -profile local
``` 

## Step 4: Run "se_resistome" and point to non host reads with `--reads`

Parameters that have to change:
* `--pipeline` ==> `--pipeline se_resistome`
* `--reads`  ==> `--reads "SE_AMR++_analysis/HostRemoval/NonHostFastq/*.fastq.gz"`

SNP confirmation and alignment deduplication is performed by default.

Example command:
```
nextflow run main_AMR++.nf --pipeline se_resistome --output SE_AMR++_analysis --reads "SE_AMR++_analysis/HostRemoval/NonHostFastq/*.fastq.gz" -profile local
``` 


## Optional - Step 5: Run "se_kraken" and point to non host reads with `--reads`

Parameters that have to change:
* `--pipeline` ==> `--pipeline se_kraken`
* `--kraken_db` ==> `--kraken_db /path/to/your/kraken_db`


Parameter still pointing to nonhost merged reads
* `--reads`  ==> `--reads 'SE_AMR++_analysis/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz'`

Example command:
```
nextflow run main_AMR++.nf --pipeline se_kraken --output SE_AMR++_analysis --reads "SE_AMR++_analysis/HostRemoval/NonHostFastq/*.fastq.gz" --kraken_db "/path/to/your/kraken_db" -profile local
``` 