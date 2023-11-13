#!/bin/bash

# Job script
echo "Job started at: `date`"
echo "Job ID: $SLURM_JOBID"

## Data sources and target dirs
chrom=${1}
target_dir=${2}
output_prefix="${target_dir}/chr${chrom}.high_corr"
bimfile_prefix=${3}
rsids_path="${4}/rsids.${chrom}.txt"

# find highly correlated variants for each snp
$PLINK \
    --bfile $bfile_prefix \
    --show-tags $rsids \
    --list-all \
    --tag-kb 3000 \
    --tag-r2 0.8 \
    --out $output_prefix

    