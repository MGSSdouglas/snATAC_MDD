import os
import pdb
import pickle

import os.path as osp
import pandas as pd
import numpy as np

from Bio.SeqUtils import GC


"""
This script contains the configuration and data preparation script for the
compute job for generating cluster-specific negative examples 
not overlapping with OCRs of the cluster
"""

## Configuration

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_ENV_DIR = os.environ.get("GKMSVM_ENV_DIR")
GKMSVM_SCRIPTS_DIR = os.environ.get("GKMSVM_SCRIPTS_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")


# source and target filepaths
tmp_dir = osp.join(GKMSVM_TMP_DIR, "5_2_candidate_MDD_SNPs")
os.makedirs(tmp_dir, exist_ok=True)

## Data preparation 

# load all lead snps and combine them
chrom2snps = {chrom:[] for chrom in range(1,22+1)}

# als et al
filepath = osp.join(prepared_data_dir, "als_et_al.gws.tsv")
snps = pd.read_csv(filepath, sep="\t")
for chrom, chrom_snps in snps[["Chromosome", "Name"]].groupby("Chromosome"):
    chrom2snps[chrom] += chrom_snps["Name"].tolist()

# levey et al
filepath = osp.join(prepared_data_dir, "levey_et_al.gws.tsv")
snps = pd.read_csv(filepath, sep="\t")
for chrom, chrom_snps in snps[["Chromosome", "Name"]].groupby("Chromosome"):
    chrom2snps[chrom] += chrom_snps["Name"].tolist()

# howard et al
filepath = osp.join(prepared_data_dir, "levey_et_al.gws.tsv")
snps = pd.read_csv(filepath, sep="\t")
for chrom, chrom_snps in snps[["Chromosome", "Name"]].groupby("Chromosome"):
    chrom2snps[chrom] += chrom_snps["Name"].tolist()

chrom2snps = {chrom[3:]:set(val) for chrom,val in chrom2snps.iteritems()}

# save lead snp rsids per chromosome
rsid_filepath_base = osp.join(tmp_dir, "lead_snps_per_chrom")
os.makedirs(rsid_filepath_base, exist_ok=True)
for chrom in chrom2snps:
    filepath = osp.join(rsid_filepath_base, f"rsids.{chrom}.txt")
    with open(filepath, "w") as f:
        for rsid in chrom2snps[chrom]:
            f.write(f"{rsid}\n")

# 1kG EUR plinkfiles
data_id = "1000G_EUR_Phase3_plinkfiles"
bimfile_dir = osp.join(GKMSVM_RAW_DATA_DIR, data_id)
bimfile_template = "1000G.EUR.QC.{}.bim"

# target filepaths
out_filepath = osp.join(tmp_dir, "gws.ld_expansion")
os.makedirs(out_filepath, exist_ok=True)

# job params
jobs = []
for chrom in chrom2snps:
    jobs.append({"Name": f"{chrom}",
                "chrom": chrom,
                "target_dir": out_filepath,
                "bimfile_path": osp.join(bimfile_dir, bimfile_template.format(chrom))
                "rsids_dir":rsid_filepath_base})


