import pdb
import os 
import sys 
import pickle
import os.path as osp

import pandas as pd 
import numpy as np 
import pyranges as pr

from tqdm import tqdm
from pprint import pprint

"""
Utility functions
"""
def aggregate_v1(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    if len(np.unique(x)) != 1:
        return ",".join([str(y) for y in x])
    else:
        return [str(y) for y in x][0]

def aggregate_v2(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    return ",".join([str(y) for y in x])

"""
Configuration
"""

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
output_plot_dir = f"{output_dir}/plots"

# load most recent pipeline_data
with open(f"{temp_data_dir}/merged_data.v11.pkl", "rb") as f:
    pipeline_data = pickle.load(f)

# iterate over cell types in pipeline_data and concatenate snp tables
tables = []
for cell_type in pipeline_data:
    snp_table = pipeline_data[cell_type]["original_snp_tables"]
    tables.append(snp_table)
snp_table = pd.concat(tables)

# retain the first row per SNP in Name column 
snp_table_ = snp_table.groupby("Name").first().reset_index()

# retain Chromosome, Start, End and Name columns 
snp_table_ = snp_table_[["Chromosome", "Start", "End", "Name"]]

# convert 0 based coordinates to 1 based
snp_table_["Start"] = snp_table_["Start"] + 1

# save the table
snp_table_.to_csv(f"{temp_data_dir}/candSNP_locs.tsv", sep="\t", index=False)

print("Script finished")
