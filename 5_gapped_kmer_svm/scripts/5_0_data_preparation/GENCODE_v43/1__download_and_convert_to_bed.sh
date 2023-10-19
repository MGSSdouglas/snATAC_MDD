#!/bin/bash

# This script downloads GENCODE_v43 to raw_data directory, and converts it to bed file
module load bedtools

# config
data_id="GENCODE_v43"
raw_data_dir="$GKMSVM_RAW_DATA_DIR/$data_id"
mkdir -p $raw_data_dir
cd $raw_data_dir

# download GENCODE_v43
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gff3.gz

# create prepared data dir
prepared_data_dir="$GKMSVM_PREPARED_DATA_DIR/$data_id"

### Convert GENCODE gff3 file to bed file

# file paths
GENCODE_gff_path="${raw_data_path}/gencode.v43.annotation.gff3.gz"
GENCODE_genePred_path="${prepared_data_dir}/gencode.v43.annotation.genePred"
GENCODE_bed_path="${prepared_data_dir}/gencode.v43.annotation.bed"
sorted_GENCODE_bed_path="${prepared_data_dir}/gencode.v43.annotation.sorted.bed"

# tool paths
ucsc_utils_basepath="$GKMSVM_BIN_DIR/ucsc_utils"
gff3ToGenePred="${ucsc_utils_basepath}/gff3ToGenePred"
genePredToBed="${ucsc_utils_basepath}/genePredToBed"

# step 1: convert gff3 to genePred
${gff3ToGenePred} ${GENCODE_gff_path} ${GENCODE_genePred_path}

# step 2: convert genePred to bed
${genePredToBed} ${GENCODE_genePred_path} ${GENCODE_bed_path}

# step3: sort GENCODE bed file
sort -k1,1 -k2,2n ${GENCODE_bed_path} > ${sorted_GENCODE_bed_path}


echo "Script finished"
