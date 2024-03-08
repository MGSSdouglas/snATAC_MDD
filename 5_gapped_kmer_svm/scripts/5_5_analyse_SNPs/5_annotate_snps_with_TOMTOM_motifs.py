import pdb
import os 
import sys 
import pickle
import os.path as osp

import pandas as pd 
import numpy as np 
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import statsmodels as sm
import pyranges as pr

from Bio import SeqIO
from tqdm import tqdm
from copy import deepcopy
from pprint import pprint
from statsmodels.graphics.gofplots import qqplot_2samples
from statsmodels.distributions.empirical_distribution import ECDF


"""
Configuration
"""

sub2broad = {}
for k, v in broad2sub.items():
    for x in v:
        sub2broad[x] = k

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"
os.makedirs(temp_data_dir, exist_ok=True)

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
output_plot_dir = f"{output_dir}/plots"
os.makedirs(output_data_dir, exist_ok=True)
null_sequences_count = 10

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

# read latest pipeline_data structure 
with open(f"{temp_data_dir}/merged_data.v6.pkl", "rb") as f:
    pipeline_data = pickle.load(f)

### Join TFBS match results to the snp table
for cell_type in pipeline_data:
    # get original snp table
    snp_table = pipeline_data[cell_type]['original_snp_tables']
    snp_table.set_index("Name", inplace=True)

    # CISBP
    ### prepare analysis cell type specific tf motif matches table
    motif_table = f"{temp_data_dir}/TFBS/tomtom_matches/{cell_type}__CISBP_unnested.tsv"
    motif_table = pd.read_csv(motif_table, sep="\t")
    motif_table = motif_table[motif_table["analysis_cell_type"] == cell_type]
    motif_table = motif_table[motif_table["match_qval"] < 0.1]
    motif_table = motif_table.groupby(["name"]).agg(aggregate_v2)
    motif_table.drop(columns=["alphabet", "strand", "icscore", "nsites", \
        "bkgsites", "pval", "qval", "eval"], inplace=True)
    motif_table.columns = ["CISBP__" + x for x in motif_table.columns]

    ### join motif matches to the snp table
    snp_table = snp_table.join(motif_table, how="left")

    # add snp table to pipeline_data
    pipeline_data[cell_type]['original_snp_tables'] = snp_table

# save the resulting pipeline_data structure 
with open(f"{temp_data_dir}/merged_data.v7.pkl", "wb") as f:
    pickle.dump(pipeline_data, f)

print("Script finished")
