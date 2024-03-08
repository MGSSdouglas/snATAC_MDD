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
Configuration
"""

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
output_plot_dir = f"{output_dir}/plots"

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

def aggregate_v3(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    return sum([y for y in x])


""" Candidate SNPs """
# read candidate snp table
candSNP_table = pd.read_csv(f"{output_data_dir}/cdSNP_table.tsv", sep="\t")

##### filling the following table 
snp_stat_table = {
    "snp__Name": [],
    "pipeline__cell_type": [],
    "pipeline__ssSNP": [],
    "peak__cell_type": [],
    "peak__is_broad_marker": [],
    "peak__is_broad_diff": [],
    "peak__is_subcluster_marker": [],
    "peak__is_subcluster_diff": [],
}
# group candidate snp table with respect to snp__Name column
groups = candSNP_table.groupby("snp__Name")
for rsid, table in groups:
    for idx, row in table.iterrows():
        
        # get peak details
        peak__origin_cell_type, peak__origin_peak_type = row['peak__origin_cell_type'], row['peak__origin_peak_type']
        if peak__origin_cell_type in broad_clusters:
            # set some details
            snp_stat_table["snp__Name"].append(rsid)
            snp_stat_table["pipeline__cell_type"].append(row["pipeline__cell_type"])
            snp_stat_table["pipeline__ssSNP"].append(row["pipeline__ssSNP"])
            snp_stat_table["peak__cell_type"].append(row["pipeline__cell_type"])
            if peak__origin_peak_type == "marker":
                snp_stat_table["peak__is_broad_marker"].append(True)
                snp_stat_table["peak__is_broad_diff"].append(False)
                snp_stat_table["peak__is_subcluster_marker"].append(False)
                snp_stat_table["peak__is_subcluster_diff"].append(False)
            elif peak__origin_peak_type == "diff":
                snp_stat_table["peak__is_broad_marker"].append(False)
                snp_stat_table["peak__is_broad_diff"].append(True)
                snp_stat_table["peak__is_subcluster_marker"].append(False)
                snp_stat_table["peak__is_subcluster_diff"].append(False)
            else:
                raise ValueError("peak__origin_peak_type {} not recognized".format(peak__origin_peak_type))
        elif peak__origin_cell_type in cell_types:
            # set some details
            snp_stat_table["snp__Name"].append(rsid)
            snp_stat_table["pipeline__cell_type"].append(row["pipeline__cell_type"])
            snp_stat_table["pipeline__ssSNP"].append(row["pipeline__ssSNP"])
            snp_stat_table["peak__cell_type"].append(peak__origin_cell_type)
            if peak__origin_peak_type == "marker":
                snp_stat_table["peak__is_broad_marker"].append(False)
                snp_stat_table["peak__is_broad_diff"].append(False)
                snp_stat_table["peak__is_subcluster_marker"].append(True)
                snp_stat_table["peak__is_subcluster_diff"].append(False)
            elif peak__origin_peak_type == "diff":
                snp_stat_table["peak__is_broad_marker"].append(False)
                snp_stat_table["peak__is_broad_diff"].append(False)
                snp_stat_table["peak__is_subcluster_marker"].append(False)
                snp_stat_table["peak__is_subcluster_diff"].append(True)
            else:
                raise ValueError("peak__origin_peak_type {} not recognized".format(peak__origin_peak_type))
        else:
            raise ValueError("peak__origin_cell_type {} not recognized".format(peak__origin_cell_type))
        
# convert to pandas dataframe
snp_stat_table = pd.DataFrame(snp_stat_table)
# create peak__is_marker column 
snp_stat_table["peak__is_marker"] = snp_stat_table["peak__is_broad_marker"] | snp_stat_table["peak__is_subcluster_marker"]
# create peak__is_diff column
snp_stat_table["peak__is_diff"] = snp_stat_table["peak__is_broad_diff"] | snp_stat_table["peak__is_subcluster_diff"]
# retain rows where pipeline cell type and peak cell type are same 
snp_stat_table = snp_stat_table[snp_stat_table["pipeline__cell_type"] == snp_stat_table["peak__cell_type"]]
# save as temp.tsv to temp data dir
snp_stat_table.to_csv(f"{output_data_dir}/candSNP.tsv", sep="\t", index=False)

##### find candidate SNPs that are exclusive to each subcluster 
cell_type_exclusive_snps = {cell_type:[] for cell_type in cell_types}
unique_snps = snp_stat_table.snp__Name.unique().tolist()
for rsid in unique_snps:
    # get table rows that are associated with this rsid 
    table = snp_stat_table[snp_stat_table["snp__Name"] == rsid]
    if table.pipeline__cell_type.unique().shape[0] == 1: # this snp is considered for the pipeline of only one cell type
        cell_type_exclusive_snps[table.pipeline__cell_type.unique()[0]].append(rsid)
    else:
        # this snp is considered for pipeline of more than one cell type
        continue

# pdb.set_trace()
pprint({key:len(value) for key,value in cell_type_exclusive_snps.items()})

# create a table from cell type exclusive snps 
cell_type_exclusive_snps_table = {
    "pipeline__cell_type": [],
    "rsids": [],
    "count": []
}
for key, value in cell_type_exclusive_snps.items():
    cell_type_exclusive_snps_table["pipeline__cell_type"].append(key)
    cell_type_exclusive_snps_table["rsids"].append(",".join(value))
    cell_type_exclusive_snps_table["count"].append(len(value))
cell_type_exclusive_snps_table = pd.DataFrame(cell_type_exclusive_snps_table)
# save this table 
cell_type_exclusive_snps_table.to_csv(f"{output_data_dir}/candSNP.cell_type_exclusive_snps.tsv", sep="\t", index=False)

##### generate SNP to eQTL mapping table 

##### process this table create cell type statistics table 
cell_type_statistics_table = {
    "cell_type": [],
    "cand_snp_count": [],
    "cand_snp_in_broad_diff":[],
    "cand_snp_in_broad_marker":[],
    "cand_snp_in_subcluster_diff":[],
    "cand_snp_in_subcluster_marker":[],
    "cand_snp_in_marker":[],
    "cand_snp_in_diff":[]
}

# add stats across all cell types 
cell_type_statistics_table["cell_type"].append("across_cell_types")
cell_type_statistics_table["cand_snp_count"].append(len(snp_stat_table["snp__Name"].unique().tolist()))
cell_type_statistics_table["cand_snp_in_broad_diff"].append(len(snp_stat_table[snp_stat_table["peak__is_broad_diff"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["cand_snp_in_broad_marker"].append(len(snp_stat_table[snp_stat_table["peak__is_broad_marker"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["cand_snp_in_subcluster_diff"].append(len(snp_stat_table[snp_stat_table["peak__is_subcluster_diff"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["cand_snp_in_subcluster_marker"].append(len(snp_stat_table[snp_stat_table["peak__is_subcluster_marker"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["cand_snp_in_marker"].append(len(snp_stat_table[snp_stat_table["peak__is_marker"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["cand_snp_in_diff"].append(len(snp_stat_table[snp_stat_table["peak__is_diff"]]["snp__Name"].unique().tolist()))

for cell_type in cell_types:
    table = snp_stat_table[snp_stat_table["pipeline__cell_type"] == cell_type]
    if len(table) == 0:
        # add empty row
        cell_type_statistics_table["cell_type"].append(cell_type)
        cell_type_statistics_table["cand_snp_count"].append(0)
        cell_type_statistics_table["cand_snp_in_broad_diff"].append(0)
        cell_type_statistics_table["cand_snp_in_broad_marker"].append(0)
        cell_type_statistics_table["cand_snp_in_subcluster_diff"].append(0)
        cell_type_statistics_table["cand_snp_in_subcluster_marker"].append(0)
        cell_type_statistics_table["cand_snp_in_marker"].append(0)
        cell_type_statistics_table["cand_snp_in_diff"].append(0)
    else:
        # add row with values
        cell_type_statistics_table["cell_type"].append(cell_type)
        # drop cell type column
        table = table.drop(columns=["pipeline__cell_type", 'pipeline__ssSNP', 'peak__cell_type'])
        # group by snps and use aggregate.v3 function 
        table = table.groupby("snp__Name").agg(aggregate_v3)

        cell_type_statistics_table["cand_snp_count"].append(len(table))
        cell_type_statistics_table["cand_snp_in_broad_diff"].append(table["peak__is_broad_diff"].sum())
        cell_type_statistics_table["cand_snp_in_broad_marker"].append(table["peak__is_broad_marker"].sum())
        cell_type_statistics_table["cand_snp_in_subcluster_diff"].append(table["peak__is_subcluster_diff"].sum())
        cell_type_statistics_table["cand_snp_in_subcluster_marker"].append(table["peak__is_subcluster_marker"].sum())
        cell_type_statistics_table["cand_snp_in_marker"].append(table["peak__is_marker"].sum())
        cell_type_statistics_table["cand_snp_in_diff"].append(table["peak__is_diff"].sum())

cell_type_statistics_table = pd.DataFrame(cell_type_statistics_table)

# save table
cell_type_statistics_table.to_csv(f"{output_data_dir}/candSNP.cell_type_statistics.tsv", sep="\t", index=False)


""" ssSNPs """

##### retain rows where pipeline ssSNP is true 
snp_stat_table = snp_stat_table[snp_stat_table["pipeline__ssSNP"] == True]

snp_stat_table.to_csv(f"{output_data_dir}/ssSNP.tsv", sep="\t", index=False)

##### find ssSNPs that are exclusive to each subcluster 
cell_type_exclusive_snps = {cell_type:[] for cell_type in cell_types}
unique_snps = snp_stat_table.snp__Name.unique().tolist()
for rsid in unique_snps:
    # get table rows that are associated with this rsid 
    table = snp_stat_table[snp_stat_table["snp__Name"] == rsid]
    if table.pipeline__cell_type.unique().shape[0] == 1: # this snp is found as statistically significant for the pipeline of only one cell type
        cell_type_exclusive_snps[table.pipeline__cell_type.unique()[0]].append(rsid)
    else:
        # this snp is found as statistically significant for pipeline of more than one cell type
        continue

# pdb.set_trace()
pprint({key:len(value) for key,value in cell_type_exclusive_snps.items()})

# create a table from cell type exclusive snps 
cell_type_exclusive_snps_table = {
    "pipeline__cell_type": [],
    "rsids": [],
    "count": []
}
for key, value in cell_type_exclusive_snps.items():
    cell_type_exclusive_snps_table["pipeline__cell_type"].append(key)
    cell_type_exclusive_snps_table["rsids"].append(",".join(value))
    cell_type_exclusive_snps_table["count"].append(len(value))
cell_type_exclusive_snps_table = pd.DataFrame(cell_type_exclusive_snps_table)
# save this table 
cell_type_exclusive_snps_table.to_csv(f"{output_data_dir}/ssSNP.cell_type_exclusive_snps.tsv", sep="\t", index=False)

##### process this table create cell type statistics table 
cell_type_statistics_table = {
    "cell_type": [],
    "ss_snp_count": [],
    "ss_snp_in_broad_diff":[],
    "ss_snp_in_broad_marker":[],
    "ss_snp_in_subcluster_diff":[],
    "ss_snp_in_subcluster_marker":[],
    "ss_snp_in_marker":[],
    "ss_snp_in_diff":[]
}

# add stats across all cell types 
cell_type_statistics_table["cell_type"].append("across_cell_types")
cell_type_statistics_table["ss_snp_count"].append(len(snp_stat_table["snp__Name"].unique().tolist()))
cell_type_statistics_table["ss_snp_in_broad_diff"].append(len(snp_stat_table[snp_stat_table["peak__is_broad_diff"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["ss_snp_in_broad_marker"].append(len(snp_stat_table[snp_stat_table["peak__is_broad_marker"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["ss_snp_in_subcluster_diff"].append(len(snp_stat_table[snp_stat_table["peak__is_subcluster_diff"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["ss_snp_in_subcluster_marker"].append(len(snp_stat_table[snp_stat_table["peak__is_subcluster_marker"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["ss_snp_in_marker"].append(len(snp_stat_table[snp_stat_table["peak__is_marker"]]["snp__Name"].unique().tolist()))
cell_type_statistics_table["ss_snp_in_diff"].append(len(snp_stat_table[snp_stat_table["peak__is_diff"]]["snp__Name"].unique().tolist()))

for cell_type in cell_types:
    table = snp_stat_table[snp_stat_table["pipeline__cell_type"] == cell_type]
    if len(table) == 0:
        # add empty row
        cell_type_statistics_table["cell_type"].append(cell_type)
        cell_type_statistics_table["ss_snp_count"].append(0)
        cell_type_statistics_table["ss_snp_in_broad_diff"].append(0)
        cell_type_statistics_table["ss_snp_in_broad_marker"].append(0)
        cell_type_statistics_table["ss_snp_in_subcluster_diff"].append(0)
        cell_type_statistics_table["ss_snp_in_subcluster_marker"].append(0)
        cell_type_statistics_table["ss_snp_in_marker"].append(0)
        cell_type_statistics_table["ss_snp_in_diff"].append(0)
    else:
        # add row with values
        cell_type_statistics_table["cell_type"].append(cell_type)
        # drop cell type column
        table = table.drop(columns=["pipeline__cell_type", 'pipeline__ssSNP', 'peak__cell_type'])
        # group by snps and use aggregate.v3 function 
        table = table.groupby("snp__Name").agg(aggregate_v3)

        cell_type_statistics_table["ss_snp_count"].append(len(table))
        cell_type_statistics_table["ss_snp_in_broad_diff"].append(table["peak__is_broad_diff"].sum())
        cell_type_statistics_table["ss_snp_in_broad_marker"].append(table["peak__is_broad_marker"].sum())
        cell_type_statistics_table["ss_snp_in_subcluster_diff"].append(table["peak__is_subcluster_diff"].sum())
        cell_type_statistics_table["ss_snp_in_subcluster_marker"].append(table["peak__is_subcluster_marker"].sum())
        cell_type_statistics_table["ss_snp_in_marker"].append(table["peak__is_marker"].sum())
        cell_type_statistics_table["ss_snp_in_diff"].append(table["peak__is_diff"].sum())

cell_type_statistics_table = pd.DataFrame(cell_type_statistics_table)

# save table
cell_type_statistics_table.to_csv(f"{output_data_dir}/ssSNP.cell_type_statistics.tsv", sep="\t", index=False)


print("Script finished")
