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

You can run the full AMR++ pipeline (without kraken) by specifying the path to your `--reads` and using the `--pipeline merged_AMR` flag. This could work well if you re-direct the results or the work directory `-w /path/to/shared_drive`.

## Step 1 : Run the "eval_qc" pipeline as normal with `--reads`

Example command:
```
nextflow run main_AMR++.nf --pipeline eval_qc --output Merged_AMR++_analysis --reads "data/raw/*_R{1,2}.fastq.gz" -profile local
``` 

## Step 2: Run "trim_qc" pipeline as normal with `--reads`
Modify trimming parameters as needed in the `params.txt` file. 


Example command:
```
nextflow run main_AMR++.nf --pipeline trim_qc --output Merged_AMR++_analysis --reads "data/raw/*_R{1,2}.fastq.gz" -profile local
``` 

## Step 3: Run new "merge_reads" pipeline as normal with `--reads`

### Parameters that have to change:
* `--reads` ==> `--reads "Merged_AMR++_analysis/QC_trimming/Paired/*{1,2}P.fastq.gz"`
* `--pipeline` ==> `--pipeline merge_reads`

This will output two files per sample, the "extendedFrags" (merged) and "notCombined" (unmerged).
Example command:
```
nextflow run main_AMR++.nf --pipeline merge_reads --output Merged_AMR++_analysis --reads "Merged_AMR++_analysis/QC_trimming/Paired/*{1,2}P.fastq.gz" -profile local
``` 


## Step 4: Run "merged_rm_host" with a new flag, `--merged_reads`

### Parameters to change

Parameters that have to change:
* `--pipeline` ==> `--pipeline merged_rm_host`
* `--merged_reads`  ==> `--merged_reads 'Merged_AMR++_analysis/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz'`
    * Remember, it's very important to use the single quote (') and not the backtick or backquote (`) that's on the same key as the tilde. 
* `host` ==> `--host "/path/to/your/host/chr21.fasta.gz"` 
    * remember, you can change this in `params.config` file or add it to your nextflow command.
    * On grace, bovine: `/scratch/group/big_scratch/SHARED_resources/host_genome/GCF_002263795.3_ARS-UCD2.0_genomic.fna`

Example command:
```
nextflow run main_AMR++.nf --pipeline merged_rm_host --output Merged_AMR++_analysis --merged_reads 'Merged_AMR++_analysis/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz' -profile local
``` 

## Step 5: Run "merged_resistome" and point to non host reads with `--merged_reads`

Parameters that have to change:
* `--pipeline` ==> `--pipeline merged_resistome`
* `--merged_reads`  ==> `--merged_reads 'Merged_AMR++_analysis/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz'`

SNP confirmation and alignment deduplication is performed by default.

Example command:
```
nextflow run main_AMR++.nf --pipeline merged_resistome --output Merged_AMR++_analysis --merged_reads 'Merged_AMR++_analysis/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz' -profile local
``` 


## Optional - Step 6: Run "merged_kraken" and point to non host reads with `--merged_reads`

Parameters that have to change:
* `--pipeline` ==> `--pipeline merged_kraken`
* `--kraken_db` ==> `--kraken_db /path/to/your/kraken_db`


Parameter still pointing to nonhost merged reads
* `--merged_reads`  ==> `--merged_reads 'Merged_AMR++_analysis/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz'`

Example command:
```
nextflow run main_AMR++.nf --pipeline merged_kraken --output Merged_AMR++_analysis --merged_reads 'Merged_AMR++_analysis/HostRemoval/NonHostFastq/*.{merged,unmerged}.non.host.fastq.gz' --kraken_db "/path/to/your/kraken_db" -profile local
``` 