Usage
-----

### Display Help Message

The `help` parameter displays the available options and commands.
```
$ nextflow run auir.nf --help
```

### File Inputs

#### Set custom sequence data

The `reads` parameter accepts sequence files in standard fastq and gz format.
```
$ nextflow run main.nf --reads "data/raw/*_R{1,2}.fastq"
```

#### Set host genome

The `host` parameter accepts a fasta formatted host genome.
```
$ nextflow run main.nf --host "data/host/bovine.fa"
```

#### Set host index

The `index` parameter allows you to upload pre-built host indexes produced by BWA.
```
$ nextflow run main.nf --host "data/host/bovine.fa" --index "data/index/*"
```

#### Set resistance database

The `amr` parameter accepts a fasta formatted resistance database. 
```
$ nextflow run main.nf --amr "data/amr/megares_database_v1.01.fasta"
```

#### Set annotation database

The `annotation` parameter accepts a csv formatted annotation database.
```
$ nextflow run main.nf --annotation "data/amr/megares_annotations_v1.01.csv"
```

#### Set adapter file

The `adapters` parameter accepts a fasta formatted adapter file.
```
$ nextflow run main.nf --adapters "data/adapters/adapters.fa"
```

### File Outputs

#### Set output directory

The `output` parameter writes output files to the specified directory.
```
$ nextflow run main.nf --output "test/"
```

### Trimming Options

#### Set custom trimming parameters

```
$ nextflow run main.nf \
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
$ nextflow run main.nf \
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
$ nextflow run main.nf --threads 8
```
