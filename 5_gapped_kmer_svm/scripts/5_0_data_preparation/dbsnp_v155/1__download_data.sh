#!/bin/bash

# This script downloads dbSNP (build 155) bigBed files from UCSC genome browser to raw data directory

# config
data_id="dbsnp155"
raw_data_dir="$GKMSVM_RAW_DATA_DIR/$data_id"
mkdir -p $raw_data_dir
cd $raw_data_dir

# download dbSNP build 155 
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb
wget http://hgdownload.soe.ucsc.edu/gbdb/hgFixed/dbSnp/dbSnp155Details.tab.gz

echo "Script finished"
