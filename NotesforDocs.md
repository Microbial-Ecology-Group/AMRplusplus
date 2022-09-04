# Notes for website


- Have to add intro tutorial with demo

- Requirements
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n AMR++ snakemake git nextflow
mamba activate AMR++


mamba create -c conda-forge -c bioconda -n AMR_env git python=3 trimmomatic multiqc bwa samtools bedtools kraken2 multiqc fastqc


# Parameter seelction
    Parameters specified on the command line (--something value)

    Parameters provided using the -params-file option

    Config file specified using the -c my_config option

    The config file named nextflow.config in the current directory

    The config file named nextflow.config in the workflow project directory

    The config file $HOME/.nextflow/config

    Values defined within the pipeline script itself (e.g. main.nf)



Bootstrap: docker
From: debian:jessie-slim

#Includes trimmomatic, samtools, bwa, bedtools, vcftools, htslib,  kraken2, SNPfinder, freebayes, bbmap

%environment
    export LC_ALL=C

%post
    apt update \
    && apt install -y --no-install-recommends \
    build-essential ca-certificates sudo tcsh\
    git make automake autoconf wget gzip unzip sed\
    zlib1g-dev curl libbz2-dev locales libncurses5-dev liblzma-dev libcurl4-openssl-dev software-properties-common apt-transport-https\
    python3-pip python3-docopt python3-pytest python-dev python3-dev\
     libssl-dev fonts-texgyre \
    gcc g++ gfortran libblas-dev liblapack-dev dos2unix libstdc++6\
    r-base-core r-recommended libgl1-mesa-glx libegl1-mesa  \
    libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
    && rm -rf /var/lib/apt/lists/*

    curl -s https://get.sdkman.io | bash
    source "$HOME/.sdkman/bin/sdkman-init.sh"
    sdk install java

    wget -c https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
    sh Anaconda3-2022.05-Linux-x86_64.sh -bfp /usr/local

    # add bioconda channels
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # install bulk of bioinformatic tools using conda
    #conda create -n AmrPlusPlus_env python=3 trimmomatic multiqc bwa samtools bedtools freebayes bbmap vcftools htslib kraken2

    conda create -c conda-forge -n mamba_base mamba
    conda activate mamba_base
    mamba create -c conda-forge -c bioconda -n AMR++ git python=3 trimmomatic multiqc bwa samtools bedtools kraken2 multiqc fastqc
    mamba activate AMR++ 


    #. /usr/local/bin/activate AmrPlusPlus_env
    
    #ln -s /usr/local/envs/AmrPlusPlus_env/bin/* /usr/local/bin/
    
    #Still experimenting with how to change $PATH location. 
    echo 'export PATH=$PATH:/usr/local/envs/AmrPlusPlus_env/bin/' >> $SINGULARITY_ENVIRONMENT

    # SNPfinder
    cd /usr/local
    git clone https://github.com/cdeanj/snpfinder.git
    cd snpfinder
    make
    cp snpfinder /usr/local/bin
    cd /

    # Make sure all the tools have the right permissions to use the tools
    chmod -R 777 /usr/local/
    
%test


JAVA_HOME = "/home/enrique/.sdkman/candidates/java/current:-:-:-:-:-:-:-"




/opt/conda/bin/conda



Bootstrap: docker
From: continuumio/miniconda3

%labels
    Version v0.0.1

%help
    This is a demo container used to illustrate a def file.
%post

    #wget -c https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
    #sh Anaconda3-2022.05-Linux-x86_64.sh -bfp /usr/local

    # add bioconda channels
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # install bulk of bioinformatic tools using conda
    conda create -c conda-forge -n mamba_base mamba
    pwd
    conda init
    conda activate mamba_base
    mamba create -c conda-forge -c bioconda -n AMR++ git python=3 trimmomatic multiqc bwa samtools bedtools kraken2 multiqc fastqc
    mamba activate AMR++ 