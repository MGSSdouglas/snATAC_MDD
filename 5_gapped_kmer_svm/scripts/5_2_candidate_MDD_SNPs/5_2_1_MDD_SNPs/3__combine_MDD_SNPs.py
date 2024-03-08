import os
import pdb
import pickle

import os.path as osp
import pandas as pd
import numpy as np

from Bio.SeqUtils import GC


"""
This script loads LD expanded SNPs and maps them to MDD GWAS lead SNPs
"""

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_ENV_DIR = os.environ.get("GKMSVM_ENV_DIR")
GKMSVM_SCRIPTS_DIR = os.environ.get("GKMSVM_SCRIPTS_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# create temporary path for current analysis
tmp_dir = osp.join(GKMSVM_TMP_DIR, "5_2_candidate_MDD_SNPs")
os.makedirs(tmp_dir, exist_ok=True)

# load lead MDD GWAS SNPs
data_id = "snps/lead_MDD_GWAS_SNPs"
prepared_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id)

# howard et al.
filepath = osp.join(prepared_data_dir, "howard_et_al.gws.tsv")
howard_et_al = pd.read_csv(filepath, sep="\t")

# als et al.
filepath = osp.join(prepared_data_dir, "als_et_al.gws.tsv")
als_et_al = pd.read_csv(filepath, sep="\t")

# levey et al.
filepath = osp.join(prepared_data_dir, "levey_et_al.gws.tsv")
levey_et_al = pd.read_csv(filepath, sep="\t")

lead_snps = pd.concat([howard_et_al, als_et_al, levey_et_al], axis=0)

# load LD expanded SNPs
filepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "ld_expand.tsv")
ld_expanded_snps = pd.read_csv(filepath, sep="\t")

# load finemapped SNPs
data_id = "snps/MDD_GWAS_finemapping"
prepared_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id, sep = "/")
filepath = osp.join(prepared_data_dir, "snps.hg38_information.tsv")
finemapped_snps = pd.read_csv(filepath, index=False, sep="\t")

# merge lead, ld expanded and finemapped snps
mdd_snps = pd.concat([lead_snps, ld_expanded_snps, finemapped_snps], axis=0)
mdd_snps = mdd_snps.drop_duplicates(subset="Name")

# save mdd snps
filepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "mdd_snps.tsv")
mdd_snps.to_csv(filepath, sep="\t", index=False)

