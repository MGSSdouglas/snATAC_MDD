import pdb 
import pickle
import os 

import os.path as osp
import numpy as np 
import pandas as pd 

from tqdm import tqdm

"""
This script prepares the following peaksets:
- 1001-bp snATAC cluster marker cCREs 
- 1001-bp broad cell type marker cCREs
- 1001-bp snATAC cluster DARs
- 1001-bp broad cell type DARs
- 1001-bp snATAC cluster OCRs
"""


# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# configure important paths
data_id = "peaks"
raw_data_dir = osp.join(GKMSVM_RAW_DATA_DIR, data_id)
processed_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id)
os.makedirs(processed_data_dir, exist_ok=True)


""" 1001-bp snATAC cluster marker cCREs """

data_dir = osp.join(raw_data_dir, "cluster__marker_cCREs")
peak_type = "marker"

# load 501-bp peaks
peaks = []
for filename in os.listdir(data_dir):

    cluster = filename.split("_marker.tsv")[0]

    filepath = osp.join(data_dir, filename)
    table = pd.read_csv(filepath, sep="\t")
    table.drop(columns=['idx', 'MeanDiff', 'width'], inplace=True)

    # add origin cell type and peak type
    table["origin_cell_type"] = cluster
    table["origin_peak_type"] = peak_type

    # add peak identifier column
    # Step 1 : process start and end locations by converting peak from 1001bp to 501bp while preserving peak center
    table["Start_501bp"] = table["start"] + 250
    table["End_501bp"] = table["end"] - 250
    # Step 2 : add peak identifier column with collowing format: 
    #    <chr>-<start>-<end>
    table["peak_id"] = table["chr"] + "-" + table["Start_501bp"].astype(str) + "-" + table["End_501bp"].astype(str)

    # add an absolute peak identifier column
    # structure: <peak_id>-<origin_cell_type>-<origin_peak_type>
    table["absolute_peak_id"] = table["peak_id"] + "-" + table["origin_cell_type"] + "-" + table["origin_peak_type"]

    peaks.append(table)

# merge all marker peak tables
peaks = pd.concat(peaks, axis=0)

peaks.rename(columns={"chr":"Chromosome", "start":"Start", "end":"End"}, inplace=True)

peaks.sort_values(by="FDR", inplace=True)
peaks.drop_duplicates(subset="absolute_peak_id", keep="first", inplace=True)

peaks["cluster_id"] = peaks["origin_cell_type"]
peaks.sort_values(by="cluster_id", inplace=True)

try:
    # sanity check: absolute peak identifier column should be unique
    assert peaks["absolute_peak_id"].nunique() == peaks.shape[0]
except AssertionError:
    pdb.set_trace()

# save merged peak table
filepath = f"{processed_data_dir}/cluster__marker_cCREs.tsv"
peaks.to_csv(filepath, sep="\t", index=False)


""" 1001-bp broad cell type marker cCREs """

data_dir = osp.join(raw_data_dir, "cell_type__marker_cCREs")
peak_type = "marker"

# load 501-bp peaks
peaks = []
for filename in os.listdir(data_dir):

    cell_type = filename.split("_marker.tsv")[0]

    filepath = osp.join(data_dir, filename)
    table = pd.read_csv(filepath, sep="\t")
    table.drop(columns=['idx', 'MeanDiff', "indices"], inplace=True)
    table.rename(columns={"seqnames":"chr"}, inplace=True)

    # add origin cell type and peak type
    table["origin_cell_type"] = cluster
    table["origin_peak_type"] = peak_type
    
    # add peak identifier column
    # Step 1 : process start and end locations by converting peak from 1001bp to 501bp while preserving peak center
    table["Start_501bp"] = table["start"] + 250
    table["End_501bp"] = table["end"] - 250
    # Step 2 : add peak identifier column with collowing format: 
    #    <chr>-<start>-<end>
    table["peak_id"] = table["chr"] + "-" + table["Start_501bp"].astype(str) + "-" + table["End_501bp"].astype(str)

     # add an absolute peak identifier column
    # structure: <peak_id>-<origin_cell_type>-<origin_peak_type>
    table["absolute_peak_id"] = table["peak_id"] + "-" + table["origin_cell_type"] + "-" + table["origin_peak_type"]
    
    peaks.append(table)

peaks = pd.concat(peaks, axis=0)

peaks.rename(columns={"chr":"Chromosome", "start":"Start", "end":"End"}, inplace=True)

peaks.sort_values(by="FDR", inplace=True)
peaks.drop_duplicates(subset="absolute_peak_id", keep="first", inplace=True)

peaks["cluster_id"] = peaks["origin_cell_type"]
peaks.sort_values(by="cluster_id", inplace=True)

try:
    # sanity check: absolute peak identifier column should be unique
    assert peaks["absolute_peak_id"].nunique() == peaks.shape[0]
except AssertionError:
    pdb.set_trace()

# save merged peak table
filepath = f"{processed_data_dir}/cell_type__marker_cCREs.tsv"
peaks.to_csv(filepath, sep="\t", index=False)


""" 1001-bp snATAC cluster DARs """
peak_type = "DAR"

data_path = osp.join(raw_data_dir, "cluster__DARs.tsv")
peaks = pd.read_csv(data_path, sep="\t")

peaks.rename(columns={"chr":"Chromosome", "start":"Start", "end":"End"}, inplace=True)
peaks.rename(columns={"501bp_start":"Start_501bp", "501bp_end":"End_501bp"}, inplace=True)

peaks["origin_peak_type"] = peaks["origin_peak_type"].replace("diff", "DAR")

filepath = f"{processed_data_dir}/cluster__DARs.tsv"
peaks.to_csv(filepath, sep="\t", index=False)


""" 1001-bp broad cell type DARs """
peak_type = "DAR"

data_path = osp.join(raw_data_dir, "cell_type__DARs.tsv")
peaks = pd.read_csv(data_path, sep="\t")

peaks.rename(columns={"chr":"Chromosome", "start":"Start", "end":"End"}, inplace=True)
peaks.rename(columns={"501bp_start":"Start_501bp", "501bp_end":"End_501bp"}, inplace=True)

peaks["origin_peak_type"] = peaks["origin_peak_type"].replace("diff", "DAR")

filepath = f"{processed_data_dir}/cell_type__DARs.tsv"
peaks.to_csv(filepath, sep="\t", index=False)


""" 1001-bp snATAC cluster OCRs """
data_dir = osp.join(raw_data_dir, "cluster__OCRs")
peak_type = "OCR"

# prepare target directory
target_dir = f"{processed_data_dir}/cluster__OCRs"
os.makedirs(target_dir, exist_ok=True)

# load 501-bp peaks
peaks = {}
for filename in os.listdir(data_dir):

    cluster = filename.split("_OCR.tsv")[0]

    filepath = osp.join(data_dir, filename)
    table = pd.read_csv(filepath, sep="\t")

    # add peak identifier column
    # Step 1 : process start and end locations by converting peak from 1001bp to 501bp while preserving peak center
    table["Start_501bp"] = table["start"] + 250
    table["End_501bp"] = table["end"] - 250
    # Step 2 : add peak identifier column with collowing format: 
    #    <chr>-<start>-<end>
    table["peak_id"] = table["chr"] + "-" + table["Start_501bp"].astype(str) + "-" + table["End_501bp"].astype(str)

    # remove or rename some columns
    table.rename(columns={"chr":"Chromosome", "start":"Start", "end":"End"}, inplace=True)
    table.rename(columns={"GC":"GC_501bp", "score":"neglog_MACS2_pval"}, inplace=True)
    table.drop(columns=["width", "replicateScoreQuantile", "groupScoreQuantile", "GroupReplicate"], inplace=True)

    try:
        # sanity check: peaks reproducible at least in 2 replicates should be retained
        assert table["Reproducibility"].min() >= 2
    except AssertionError:
        pdb.set_trace()

    # add origin cell type and peak type
    table["origin_cell_type"] = cluster
    table["origin_peak_type"] = peak_type

    # add an absolute peak identifier column
    # structure: <peak_id>-<origin_cell_type>-<origin_peak_type>
    table["absolute_peak_id"] = table["peak_id"] + "-" + table["origin_cell_type"] + "-" + table["origin_peak_type"]

    # record peak table
    peaks[cluster] = table

    # transfer to prepared data directory
    filepath = osp.join(target_dir, f"{cluster}_OCRs.tsv")
    table.to_csv(filepath, sep="\t", index=False)

# save prepared data structure
filepath = f"{processed_data_dir}/cluster__OCRs.pkl"
with open(filepath, "wb") as f:
    pickle.dump(peaks, f)

# merge across clusters
peaks = pd.concat(peaks.values(), axis=0)

# sort by cluster id
peaks["cluster_id"] = peaks["origin_cell_type"]
peaks.sort_values(by="cluster_id", inplace=True)

# save merged peak table
filepath = f"{processed_data_dir}/cluster__OCRs.tsv"
peaks.to_csv(filepath, sep="\t", index=False)


print("Script finished")
