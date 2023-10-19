#!/bin/bash

# This script downloads MDD GWAS sumstats to raw data directory
# NOTE: genome-wide significant SNPs for each study are available in respective supplementary tables of the original publications
# Als et al: Table S2 in https://www.medrxiv.org/content/10.1101/2022.08.24.22279149v1.supplementary-material
# Howard et al: Supplementary Table 1 in https://www.nature.com/articles/s41593-018-0326-7#Sec30
# Levey et al: Supplementary Information in https://www.nature.com/articles/s41593-021-00860-2#Sec31

### config 

# available environment variables
# $GKMSVM_WORKSPACE_DIR
# $GKMSVM_RAW_DATA_DIR
# $GKMSVM_PREPARED_DATA_DIR
# $GKMSVM_MODEL_DIR
# $GKMSVM_TMP_DIR
# $GKMSVM_BIN_DIR


data_id="MDD_GWAS/als_et_al"
data_dir="${GKMSVM_RAW_DATA_DIR}/${data_id}"
mkdir -p $data_dir
cd $data_dir
wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz 
wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/Readme%20-%20Depression%20GWAS%20%282023%29.pdf

data_id="MDD_GWAS/howard_et_al"
data_dir="${GKMSVM_RAW_DATA_DIR}/${data_id}"
mkdir -p $data_dir
cd $data_dir
wget https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y


