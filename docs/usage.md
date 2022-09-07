Usage
-----

### Display Help Message

The `help` parameter displays the available options and commands.
```
$ nextflow run main_AmrPlusPlus_v2.nf --help
```

## Parameter selection
    Parameters specified on the command line (--something value)

    Parameters provided using the -params-file option

    Config file specified using the -c my_config option

    The config file named nextflow.config in the current directory

    The config file named nextflow.config in the workflow project directory

    The config file $HOME/.nextflow/config

    Values defined within the pipeline script itself (e.g. main.nf)




### File Inputs

#### Set custom sequence data

The `reads` parameter accepts sequence files in standard fastq and gz format.
```
$ nextflow run main_AmrPlusPlus_v2.nf --reads "data/raw/*_R{1,2}.fastq"
```

#### Set host genome

The `host` parameter accepts a fasta formatted host genome.
```
$ nextflow run main_AmrPlusPlus_v2.nf --host "data/host/chr21.fasta.gz"
```

#### Set host index

The `host_index` parameter allows you to upload pre-built host indexes produced by BWA.
```
$ nextflow run main_AmrPlusPlus_v2.nf --host "data/host/chr21.fasta.gz" --host_index "data/index/*"
```

#### Set resistance database

The `amr` parameter accepts a fasta formatted resistance database. 
```
$ nextflow run main_AmrPlusPlus_v2.nf --amr "data/amr/megares_database_v1.02.fasta"
```

#### Set annotation database

The `annotation` parameter accepts a csv formatted annotation database.
```
$ nextflow run main_AmrPlusPlus_v2.nf --annotation "data/amr/megares_annotations_v1.02.csv"
```

#### Set adapter file

The `adapters` parameter accepts a fasta formatted adapter file.
```
$ nextflow run main_AmrPlusPlus_v2.nf --adapters "data/adapters/adapters.fa"
```

### File Outputs

#### Set output and work directories

The `output` parameter writes the results to the specified directory. As a nextflow variable, the `work` parameter only requires one dash and determines where the temporary files will be directed. Upon completing the run, you can delete the temporary file directory.
```
$ nextflow run main_AmrPlusPlus_v2.nf --output "test/" -work "work_dir/"
```

### Resume a pipeline run

If the pipeline run is cancelled or stopped for whatever reason, using the same command with the addition of the `-resume` flag will attempt to pick up where the pipeline stopped. 
```
$ nextflow run main_AmrPlusPlus_v2.nf --output "test/" -work "work_dir/" -resume
```

### Trimming Options

#### Set custom trimming parameters

```
$ nextflow run main_AmrPlusPlus_v2.nf \
    --reads "data/raw/*_R{1,2}.fastq" \
    --leading 3 \
    --trailing 3 \
    --minlen 36 \
    --slidingwindow 4 \
    --adapters "data/adapters/nextera.fa"
    --output "test/"
```

### Algorithm Options

#### Set custom algorithm options

```
$ nextflow run main_AmrPlusPlus_v2.nf \
    --reads "data/raw/*_R{1,2}.fastq" \
    --threshold 80 \
    --min 1 \
    --max 100 \
    --samples 5 \
    --skip 5 \
    --output "test/"
```

#### Set number of threads to use for each process

```
$ nextflow run main_AmrPlusPlus_v2.nf --threads 8
```
