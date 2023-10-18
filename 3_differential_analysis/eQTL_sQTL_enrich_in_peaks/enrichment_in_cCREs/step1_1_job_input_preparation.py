import os 
import pdb 
import sys  
import warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use('agg')

import os.path as osp
import numpy as np
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns 

from tqdm import tqdm
from math import ceil
from scipy.stats import binomtest



'''
Configuration
'''

raw_data_dir = "./data/raw"
sumstats_dir = "./data/raw/sumstats"
peaks_dir = "./data/raw/peaks"

prepared_peaks_dir = "./data/processed/peaks"
prepared_sumstats_dir = "./data/processed/sumstats"
prepared_sumstats_allsnps_dir = "./data/processed/sumstats/all_snps"

processed_data_dir = "./data/processed"
results_dir = "./results"

plots_dir = "./plots/qtl_in_peaks"
os.makedirs(plots_dir, exist_ok=True)


'''
Prepare input parameters
'''

peakset_names = []
for peakset_name in ["ct_M_0.05", "ct_M_0.05__splitted__", "cluster_M_0.05", "cluster_M_0.05__splitted__", "ct_D_0.2", "ct_D_0.2__splitted__", "cluster_D_0.2", "cluster_D_0.2__splitted__", "cluster_D_0.05", "cluster_D_0.05__splitted__"]:

    if peakset_name == "ct_M_0.05":

        ###### cell type marker peaks 
        filepath =  osp.join(prepared_peaks_dir, "broad_marker.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

        peakset_names.append("ct_M_0.05")

    elif peakset_name.startswith("ct_M_0.05__splitted__"):

        ###### cell type marker peaks 
        filepath =  osp.join(prepared_peaks_dir, "broad_marker.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        # retain cell type of interest
        for context in table["cell_type"].unique():
            peakset_names.append(f"ct_M_0.05__splitted__{context}")
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)


    elif peakset_name == "cluster_M_0.05":

        ###### cluster marker peaks (0.05, splitted) 
        filepath =  osp.join(prepared_peaks_dir, "subcluster_marker.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

        peakset_names.append("cluster_M_0.05")

    elif peakset_name.startswith("cluster_M_0.05__splitted__"):

        ###### cluster marker peaks (0.05, splitted) 
        filepath =  osp.join(prepared_peaks_dir, "subcluster_marker.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        # retain cell type of interest
        for context in table["cell_type"].unique():
            peakset_names.append(f"cluster_M_0.05__splitted__{context}")
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)


    elif peakset_name == "ct_D_0.2":

        ###### cell type DARs (0.2) 
        filepath =  osp.join(prepared_peaks_dir, "broad_diff.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        # retain cell type of interest
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

        peakset_names.append("ct_D_0.2")

    elif peakset_name.startswith("ct_D_0.2__splitted__"):

        ###### cluster DARs (0.2) 
        filepath =  osp.join(prepared_peaks_dir, "broad_diff.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        # retain cell type of interest
        for context in table["cell_type"].unique():
            peakset_names.append(f"ct_D_0.2__splitted__{context}")
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

    elif peakset_name == "cluster_D_0.2":

        ###### cluster DARs (0.2) 
        filepath =  osp.join(prepared_peaks_dir, "subcluster_diff.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

        peakset_names.append("cluster_D_0.2")

    elif peakset_name.startswith("cluster_D_0.2__splitted__"):

        ###### cluster DARs (0.2)
        filepath =  osp.join(prepared_peaks_dir, "subcluster_diff.tsv")
        table = pd.read_csv(filepath, sep="\t")
        # peaks are 1001bp, convert to 501bp
        table["Start"] = table["Start"] + 250
        table["End"] = table["End"] - 250
        # get peakName from coords
        table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
        # retain cell type of interest
        for context in table["cell_type"].unique():
            peakset_names.append(f"cluster_D_0.2__splitted__{context}")
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

    elif peakset_name == "cluster_D_0.05":

        ###### cluster DARs (0.05) 
        
        filepath = osp.join(raw_data_dir, "latest_DARs/atac_fdr5.csv")
        table = pd.read_csv(filepath, sep=",")
        table = table[["gene", "cluster_id", "p_adj.loc_treatment"]]
        table[["Chromosome", "Start", "End"]] = \
            table["gene"].str.split("-", expand=True)
        table["Start"] = table["Start"].astype(int)
        table["End"] = table["End"].astype(int)
        table.rename(columns={"gene": "peakName", "cluster_id": "cluster"}, inplace=True)
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)

        peakset_names.append("cluster_D_0.05")

    elif peakset_name.startswith("cluster_D_0.05__splitted__"):

        ###### cluster DARs (0.05) 
      
        filepath = osp.join(raw_data_dir, "latest_DARs/atac_fdr5.csv")
        table = pd.read_csv(filepath, sep=",")
        table = table[["gene", "cluster_id", "p_adj.loc_treatment"]]
        table[["Chromosome", "Start", "End"]] = \
            table["gene"].str.split("-", expand=True)
        table["Start"] = table["Start"].astype(int)
        table["End"] = table["End"].astype(int)
        table.rename(columns={"gene": "peakName", "cluster_id": "cluster"}, inplace=True)
        # retain cell type of interest
        for context in table["cluster"].unique():
            peakset_names.append(f"cluster_D_0.05__splitted__{context}")
        peaks_df = table
        peaks_pr = pr.PyRanges(table, int64=True)
    else:
        print("Invalid peakset name")
        sys.exit()
    
# NOTE: order of the following operations are important 
qtl_names = ["eQTL"] * len(peakset_names) + ["sQTL"] * len(peakset_names)
peakset_names = peakset_names + peakset_names

"""
Create job parameters dictionary
"""

# out dir
out_dir = "./step14.output"
os.makedirs(out_dir, exist_ok=True)

# see which output files are already present 
out_filenames = os.listdir(out_dir)

out_ids = []
for out_ in out_filenames:
    out_ = ".".join(out_.split(".")[:-2])
    out_qtl_type = out_.split("__")[0]
    if "splitted" in out_:
        out_peakset_name = out_.split("__")[1] + "__splitted__" + out_.split("__")[3]
    else:
        out_peakset_name = out_.split("__")[1]
    out_id = "_".join([out_peakset_name, out_qtl_type])
    out_ids.append(out_id)

job_params = [
    {
        "Name": f"{peak}_{qtl}",
        "peakset_name": peak,
        "qtl_name": qtl,
        "out_dir": out_dir
    } 
    for peak, qtl in zip(peakset_names, qtl_names)
]

print(f"Total number of jobs: {len(job_params)}")
print(f"Total jobs done: {len(out_ids)}")

job_params = [x for x in job_params if x["Name"] not in out_ids]

print(f"Number of jobs to be submitted: {len(job_params)}")
