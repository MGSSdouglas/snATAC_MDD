#!/bin/bash

### config 

# available environment variables
# $GKMSVM_WORKSPACE_DIR
# $GKMSVM_RAW_DATA_DIR
# $GKMSVM_PREPARED_DATA_DIR
# $GKMSVM_MODEL_DIR
# $GKMSVM_TMP_DIR
# $GKMSVM_BIN_DIR

### download DLPFC eqtl and sqtl data
data_id="eqtl_catalogue"
data_dir="${GKMSVM_RAW_DATA_DIR}/${data_id}"
cd $data_dir

# create sumstats dir 
mkdir -p sumstats
cd sumstats

# BrainSeq - gene expression
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000005/QTD000051/QTD000051.all.tsv.gz
echo "BrainSeq - gene expression done"

# BrainSeq - splicing
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000005/QTD000055/QTD000055.cc.tsv.gz
echo "BrainSeq - splicing done"

# CommonMind - gene expression
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000008/QTD000075/QTD000075.all.tsv.gz
echo "CommonMind - gene expression done"

# CommonMind - splicing
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000008/QTD000079/QTD000079.cc.tsv.gz
echo "CommonMind - splicing done"

# GTEx - gene expression
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000176/QTD000176.all.tsv.gz
echo "GTEx - gene expression done"

# GTEx - splicing
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000180/QTD000180.cc.tsv.gz
echo "GTEx - splicing done"

# ROSMAP - gene expression
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz
echo "ROSMAP - gene expression done"

# ROSMAP - splicing
wget http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000438/QTD000438.cc.tsv.gz
echo "ROSMAP - splicing done"

# GTEx v8 imported data - gene expression 
wget ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Frontal_Cortex_BA9.tsv.gz
echo "GTEx v8 imported data - gene expression done"

echo "Script finished"