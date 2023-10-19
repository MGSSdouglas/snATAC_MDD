#!/bin/bash

# This script downloads CIS-BP 2.0 (from MEME suite) to raw data directory

### config 

# available environment variables
# $GKMSVM_WORKSPACE_DIR
# $GKMSVM_RAW_DATA_DIR
# $GKMSVM_PREPARED_DATA_DIR
# $GKMSVM_MODEL_DIR
# $GKMSVM_TMP_DIR
# $GKMSVM_BIN_DIR

data_id="CISBP2"
data_dir="${GKMSVM_RAW_DATA_DIR}/${data_id}"
mkdir -p $data_dir
cd $data_dir

# download and extract all motifs from meme suite
wget https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.23.tgz
tar -xzf motif_databases.12.23.tgz

# locate CISBP-2.0 Homo sapiens motifs
mv motif_databases.12.23/CIS-BP_2.00/Homo_sapiens.meme ./

# delete rest of the motifs 
rm -r motif_databases.12.23


echo "Script finished"
