# Install nextflow if you don't want to install conda
# curl -s https://get.nextflow.io | bash


# Download minikraken database and unzip
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar -xvzf minikraken_8GB_202003.tgz


# Install resistome-related software

git clone https://github.com/Isabella136/AmrPlusPlus_SNP.git

git clone https://github.com/cdeanj/rarefactionanalyzer.git rarefaction_dl
cd rarefaction_dl && make
chmod 777 rarefaction 


git clone https://github.com/cdeanj/resistomeanalyzer.git resistome_dl
cd resistome_dl && make
chmod 777 resistome 
