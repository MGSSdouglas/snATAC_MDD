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

# load candidate SNP annotation data 
fpath = osp.join(temp_data_dir, "annotatR_annotated.candSNP_locs.tsv")
cand_snp_annotations = pd.read_csv(fpath, sep="\t")

# rename annotations
rename_dct = {"hg38_genes_promoters": "Promoter",
	"hg38_genes_cds": "CDS",
	"hg38_genes_5UTRs": "5\'UTR",
	"hg38_genes_exons": "Exon",
	"hg38_genes_introns": "Intron",
	"hg38_genes_3UTRs": "3\'UTR",
	"hg38_genes_intergenic": "Intergenic"}
cand_snp_annotations["annot.type"] = cand_snp_annotations["annot.type"].map(rename_dct)

# group candidate SNP annotations by mcols column 
cand_snp_annotations = cand_snp_annotations.groupby("mcols").agg(aggregate_v2)
# drop these columns ["seqnames", "start", "end", "width", "strand"]
cand_snp_annotations = cand_snp_annotations.drop(columns=["seqnames", "start", "end", "width", "strand"])
# reset index 
cand_snp_annotations = cand_snp_annotations.reset_index()
# rename mcols column to Name
cand_snp_annotations = cand_snp_annotations.rename(columns={"mcols": "Name"})
# set Name column as index
cand_snp_annotations = cand_snp_annotations.set_index("Name")
# add "annotatr__" prefix to all columns 
cand_snp_annotations.columns = ["annotatr__" + str(col) for col in cand_snp_annotations.columns]

# iterate over cell types in pipeline_data 
for cell_type in pipeline_data:
    # load original_snp_tables
    snp_table = pipeline_data[cell_type]["original_snp_tables"]
    # set Name column as index 
    snp_table = snp_table.set_index("Name")

    # join snp table with cand_snp_annotations
    snp_table = snp_table.join(cand_snp_annotations, how="left")

    unique_snps = pipeline_data[cell_type]["original_snp_tables"].Name.tolist()
    assert snp_table.index.tolist() == unique_snps

    # update pipeline_data
    pipeline_data[cell_type]["original_snp_tables"] = snp_table

# load most recent pipeline_data
with open(f"{temp_data_dir}/merged_data.v12.pkl", "wb") as f:
    pickle.dump(pipeline_data, f)


print("Script finished")
