#!/bin/bash

# This script downloads hg38 genome to raw data directory

# config
data_id="hg38/grch38_p13"
data_dir="$GKMSVM_RAW_DATA_DIR/$data_id"
mkdir -p $data_dir
cd $data_dir

# download sequence information of chromosomes one by one 
# from here https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz

# extract the downloaded files
gunzip hg38.p13.fa.gz

echo "Script finished"
