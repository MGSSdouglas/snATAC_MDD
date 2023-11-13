import pdb 
import pickle
import os 

import os.path as osp
import numpy as np 
import pandas as pd 
import pyranges as pr

from tqdm import tqdm


"""
This script generates candidate MDD SNP set (cdSNPs) by overlapping MDD SNPs
with 1001-bp marker cCREs and DARs for each snATAC cluster and broad cell type
"""

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# load snATAC cluster to broad cell type mapping
filepath = osp.join(GKMSVM_RAW_DATA_DIR, "cluster2celltype.tsv")
mapping = pd.read_csv(filepath, sep="\t")
mapping = {clust:cell_type for clust, cell_type in zip(mapping["cluster"], mapping["cell_type"])}

cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")

'''
Load and prepare target SNPs and peaks
'''

# snps
filepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "mdd_snps.tsv")
snps = pd.read_csv(filepath, sep="\t", index=False)
snps.dropna(axis=0, inplace=True)
snps.drop_duplicates(subset=["Name"], keep="first")
snps["Start"] = snps["Start"].astype(int)
snps["End"] = snps["End"].astype(int)
snps["A1"] = snps["A1"].str.upper()
snps["A2"] = snps["A2"].str.upper()

# peaks
data_id = "peaks"
prepared_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id)

filepath = f"{prepared_data_dir}/cluster__marker_cCREs.tsv"
peaks = [pd.read_csv(filepath, sep="\t")]

filepath = f"{prepared_data_dir}/cell_type__marker_cCREs.tsv"
peaks += [pd.read_csv(filepath, sep="\t")]

filepath = f"{prepared_data_dir}/cluster__marker_DARs.tsv"
peaks += [pd.read_csv(filepath, sep="\t")]

filepath = f"{prepared_data_dir}/cell_type__marker_DARs.tsv"
peaks += [pd.read_csv(filepath, sep="\t")]

peaks = pd.concat(peaks, axis=0)

# split cluster and cell type peaks 

cluster_peaks = peaks[peaks["origin_cell_type"].str.len() > 3]
cell_type_peaks = peaks[peaks["origin_cell_type"].str.len() <= 3]

'''
overlap snps with peaks
'''

pr_snps = pr.PyRanges(snps, int64=True)
pr_cluster_peaks = pr.PyRanges(cluster_peaks, int64=True)
pr_cell_type_peaks = pr.PyRanges(cell_type_peaks, int64=True)

# overlap snps with cell type peaks 
snps_in_cell_type_peaks = pr_snps.join(pr_cell_type_peaks, how="left").df
snps_in_cell_type_peaks = snps_in_cell_type_peaks[snps_in_cell_type_peaks["origin_cell_type"] != "-1"]

# overlap snps with cluster peaks
snps_in_cluster_peaks = pr_snps.join(pr_cluster_peaks, how="left").df
snps_in_cluster_peaks = snps_in_cluster_peaks[snps_in_cluster_peaks["origin_cell_type"] != "-1"]

# for each cluster, combine snps in cell type and cluster peaks, and save 
snps_per_cluster = {}
for cluster in snps_in_cluster_peaks["origin_cell_type"].unique():

    cell_type = mapping[cluster]
    temp = [snps_in_cluster_peaks[snps_in_cluster_peaks["origin_cell_type"] == cluster]]
    temp += [cell_type_peaks[cell_type_peaks["origin_cell_type"] == cell_type]]
    temp = pd.concat(temp, axis=0)

    filepath = osp.join(cdsnp_basepath, cluster)
    os.makedirs(filepath, exist_ok=True)
    filepath = osp.join(filepath, "cdsnps.all.tsv")
    temp.to_csv(filepath, sep="\t", index=False)

    # retain first occurence of each snp
    temp = temp.drop_duplicates(subset=["Name"], keep="first")
    filepath = osp.join(cdsnp_basepath, cluster, "cdsnps.unique.tsv")
    temp.to_csv(filepath, sep="\t", index=False)


