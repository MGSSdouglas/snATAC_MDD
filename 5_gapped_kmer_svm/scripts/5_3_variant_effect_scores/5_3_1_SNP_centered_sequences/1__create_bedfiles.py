import pdb 
import pickle
import os 

import os.path as osp
import numpy as np 
import pandas as pd 
import pyranges as pr

from tqdm import tqdm


"""
This script generates 51bp and 201bp sequences centered at cdSNPs for each snATAC cluster
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

# this is the folder which contains cluster specific subdirectories consisting of cdSNPs
cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")

# get clusters having cdSNPs
clusters = os.listdir(cdsnp_basepath)

# iterate over clusters 
for cluster in clusters:

    # load snps for the peaks
    try:
        snps_in_peaks = pd.read_csv( \
            osp.join(*[cdsnp_basepath, cluster, "cdsnps.unique.tsv"]), sep="\t"
        )
    except pd.errors.EmptyDataError:
        print(f"EmptyDataError: {cluster}")
        continue


    '''
    Extract 51bp surrounding snps and convert the coordinates to bed format
    '''
    seqlen = 51
    offset = seqlen // 2
    fpath = osp.join(*[cdsnp_basepath, cluster,f"{seqlen}bp_seqs"])
    os.makedirs(fpath, exist_ok=True)

    table = snps_in_peaks.copy()
    table["Start"] = table["Start"] - offset
    table["End"] = table["End"] + offset
    table = table[["Chromosome", "Start", "End", "Name"]]
    pr_table = pr.PyRanges(table, int64=True)
    pr_table.to_bed(f"{fpath}/cdsnps.unique.{seqlen}bp.bed", keep = False)


    '''
    Extract 201bp surrounding snps and convert the coordinates to bed format
    '''
    seqlen = 201
    offset = seqlen // 2
    fpath = osp.join(*[cdsnp_basepath, cluster,f"{seqlen}bp_seqs"])
    os.makedirs(fpath, exist_ok=True)

    table = snps_in_peaks.copy()
    table["Start"] = table["Start"] - offset
    table["End"] = table["End"] + offset
    table = table[["Chromosome", "Start", "End", "Name"]]
    pr_table = pr.PyRanges(table)
    pr_table.to_bed(f"{fpath}/cdsnps.unique.{seqlen}bp.bed", keep = False)


