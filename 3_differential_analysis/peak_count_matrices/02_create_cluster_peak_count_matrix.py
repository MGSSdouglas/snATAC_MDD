import pdb 
import os 
import os.path as osp
import anndata as ad 
from scipy.io import mmread
from scipy.sparse import csc_matrix
from math import floor
import pandas as pd 
import numpy as np
import scanpy as sc


rawdata_path = "../raw_data"
prepared_data_path = "../prepared_data"
analysis_data_path = osp.join(rawdata_path, "merged_case_control.v4")

"""
load adatas
"""

# cases
print("> Loading case data...")
filepath = osp.join(rawdata_path, "case_adata.h5ad")
case_adata = ad.read(filepath)
print("case data shape", case_adata.shape)

# controls 
print("> Loading control data...")
filepath = osp.join(rawdata_path, "control_adata.h5ad")
control_adata = ad.read(filepath)
print("control data shape", control_adata.shape)

# compare columns of obs of adatas 
case_obs_cols = set(case_adata.obs.columns.tolist())
control_obs_cols = set(control_adata.obs.columns.tolist())
print("case - control", len(case_obs_cols - control_obs_cols))
print("control - case", len(control_obs_cols - case_obs_cols))
print("case & control", len(case_obs_cols & control_obs_cols))

# concatenate adatas
print("> Concatenating adatas ...")
concatted_adata = ad.concat([case_adata, control_adata], axis=0)
print("concatted data shape", concatted_adata.shape)

del case_adata, control_adata

# remove cells which contain atac_A32B32
print("> Removing cells whose source subject is atac_A32B32female ...")
concatted_adata = concatted_adata[~concatted_adata.obs.Subject.str.contains("atac_A32B32female")]
print("concatted data (subject removed) shape", concatted_adata.shape)


print("> Filtering peaks ... (0.005)")
cell_count = concatted_adata.shape[0]
percent_thr = 0.005
count_thr = floor(cell_count * percent_thr)
sc.pp.filter_genes(concatted_adata, min_cells=count_thr)
print(f"concatted data (cell count > {count_thr} == {percent_thr*100}% cells) shape", concatted_adata.shape)
print("> Comparing with available data ...")

#save
concatted_adata.write(osp.join("/home/anjali5/projects/def-cnagy/For_Wenmin/For_Doruk/merged_case_control.v4/", “all_adata_counts.005.mtx”))

print("Script finished")

pdb.set_trace()

