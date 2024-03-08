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
from statsmodels.graphics.gofplots import qqplot_2samples
from tqdm import tqdm
from copy import deepcopy
from pprint import pprint



"""
Configuration
"""

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"
os.makedirs(temp_data_dir, exist_ok=True)

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
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

"""
Join merged data with peaks across cell types
"""

# load merged data structure 
with open(f"{temp_data_dir}/merged_data.v1.pkl", "rb") as f:
    data = pickle.load(f)


# load peaks data structure 
with open(f"{temp_data_dir}/peaks.pkl", "rb") as f:
    peaks = pickle.load(f)
 
# merge peaks across broad and subclusters
merged_peaks = []
for key in peaks:
    for cell_type in peaks[key]:
        merged_peaks.append(peaks[key][cell_type])
merged_peaks = pd.concat(merged_peaks, axis=0)
merged_peaks = pr.PyRanges(merged_peaks, int64=True)

# join original snp table with the merged peaks 
for cell_type in data:
    snp_table = deepcopy(data[cell_type]["original_snp_tables"])
    pr_snp_table = pr.PyRanges(snp_table, int64=True)
    snp_join_peak = pr_snp_table.join(merged_peaks, how="left", suffix="__peak").df
    # sanity check (1)
    # all snps should intersect with at least one peak
    assert (snp_join_peak["origin_peak_type"] == "-1").sum() == 0, "Sanity check (1) failed"
    # group joined table with respect to Name column
    snp_join_peak["groupby_col"] = snp_join_peak["Name"]
    snp_join_peak_groups = snp_join_peak.groupby("groupby_col").agg(aggregate_v2)

    # post processing
    # drop snp related columns except rsid 
    snp_join_peak_groups.drop(columns=["Start", "End", "A1", "A2", "snp_source"], inplace=True)
    # split "Name" column by , and take the first element 
    snp_join_peak_groups["Name"] = snp_join_peak_groups["Name"].apply(lambda x: x.split(",")[0])
    # rename peak related columns 
    snp_join_peak_groups.rename(columns={
        "Chromosome": "peak__Chromosome", 
        "Start__peak": "peak__Start",
        "End__peak": "peak__End", 
        "peak_id": "peak__peak_id", 
        "501bp_start": "peak__501bp_start",
        "501bp_end": "peak__501bp_end",
        "absolute_peak_id": "peak__absolute_peak_id", 
        "origin_cell_type": "peak__origin_cell_type", 
        "origin_peak_type": "peak__origin_peak_type", 
        "peak_source": "peak__peak_source"}, inplace=True)
    # sanity check (2)
    assert set(snp_table["Name"]) == set(snp_join_peak_groups["Name"]), "Sanity check (2) failed"

    # join original snp table with peak enhanced table 
    # set index to Name column 
    snp_table = snp_table.set_index("Name")
    snp_join_peak_groups.set_index("Name", inplace=True)
    temp = snp_table.join(snp_join_peak_groups, on="Name", how="left")
    temp.reset_index(inplace=True)
    snp_table.reset_index(inplace=True)
    # sanity check (3)
    assert temp["Name"].tolist() == data[cell_type]["original_snp_tables"]["Name"].tolist(), "Sanity check (3) failed"
    
    # save to data structure
    data[cell_type]["original_snp_tables"] = temp

'''
Annotate SNPs with their sources
'''

# load initial SNP table 
target_snpset_basepath = "/home/dcakma3/scratch/mdd-prepare_target_snps/processed_data"
target_snpset_names = [
    "als_et_al__index",
    "als_et_al__index__1kg_ld_expand__0.8",
    "als_et_al__index__ld_friends__0.8",
    "levey_et_al__lead",
    "levey_et_al__lead__1kg_ld_expand__0.8",
    "howard_et_al__lead",
    "howard_et_al__lead__1kg_ld_expand__0.8",
    "howard_et_al__finemap__ALL",
]
target_snpset_filenames = ["snps.hg38_information.tsv"] * len(target_snpset_names)
snp_tables = []
for name, filename in zip(target_snpset_names, target_snpset_filenames):
    snp_table = pd.read_csv(osp.join(target_snpset_basepath, name, filename), sep="\t")
    snp_table["snp_source"] = name

    # input specific pre-processing
    if name == "als_et_al__index":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "als_et_al__index__1kg_ld_expand__0.8":
        # drop column Chromosome from snp table 
        snp_table.drop(columns=["Chromosome"], inplace=True)
        snp_table.rename(columns={"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"}, inplace=True)
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "als_et_al__index__ld_friends__0.8":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "levey_et_al__lead":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "levey_et_al__lead__1kg_ld_expand__0.8":
        snp_table.drop(columns=["Chromosome"], inplace=True)
        snp_table.rename(columns={"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"}, inplace=True)
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__lead":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__lead__1kg_ld_expand__0.8":
        snp_table.drop(columns=["Chromosome"], inplace=True)
        snp_table.rename(columns={"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"}, inplace=True)
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__ALL":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__hard_threshold_0.01":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__hard_threshold_0.005":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__hard_threshold_0.001":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__10%_cred_set":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__25%_cred_set":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]
    elif name == "howard_et_al__finemap__50%_cred_set":
        snp_table = snp_table[["Name", "Chromosome", "Start", "End", "snp_source", "A1", "A2"]]

    # remove rows where one of Start or End columns is NaN
    snp_table = snp_table[~snp_table["Start"].isna()]
    snp_table = snp_table[~snp_table["End"].isna()]

    snp_tables.append(snp_table)

# create a dictionary out of the list of snp tables
snp_sources = {}
for snpset, snp_table in zip(target_snpset_names, snp_tables):
    snp_sources[snpset] = snp_table.Name.unique().tolist()

def annotate_with_snp_source(rsid, snp_sources):
    """This function annotates snps with respective sources
    """
    annotation = []
    for snpset, snps in snp_sources.items():
        if rsid in snps:
            annotation.append(snpset)
    return ",".join(annotation)

for cell_type in data:
    data[cell_type]["original_snp_tables"]["snp_source"] = \
        data[cell_type]["original_snp_tables"]["Name"].apply(
            lambda x: annotate_with_snp_source(x, snp_sources)
        )

"""
Annotate SNPs with chromatin accessibility scores 
"""
for cell_type in data:
    # sanity check (4)
    assert data[cell_type]["original_snp_tables"].shape[0] == \
        data[cell_type]["original_gkmexplain_ve_scores"].shape[0] == \
        data[cell_type]["original_deltaSVM_ve_scores"].shape[0] == \
        data[cell_type]["original_ISM_ve_scores"].shape[0], "Sanity check (4) failed"
    # record gkmexplain variant effect score 
    data[cell_type]["original_snp_tables"]["gkmexplain_ve_score"] = list(data[cell_type]["original_gkmexplain_ve_scores"])
    # record deltaSVM variant effect score
    data[cell_type]["original_snp_tables"]["deltaSVM_ve_score"] = list(data[cell_type]["original_deltaSVM_ve_scores"])
    # record ISM variant effect score
    data[cell_type]["original_snp_tables"]["ISM_ve_score"] = list(data[cell_type]["original_ISM_ve_scores"])
    
"""
Annotate SNPs with PIPs
"""
extract_sparsepro_results = False
if extract_sparsepro_results:
    # data path
    raw_data_filepath_base = f"/home/dcakma3/scratch/mdd-prepare_target_snps/raw_data/howard_et_al__finemap"
    finemap_base_filepath = f"{raw_data_filepath_base}/MDD"
    rsid_mapping_basepath = f"{raw_data_filepath_base}/MDD/ukb"
    # pip and effect size paths
    pip_and_effect_size_path = f"{finemap_base_filepath}/all.pip"
    pips_and_effect_sizes = pd.read_csv(pip_and_effect_size_path, header=None, sep="\t")
    pips_and_effect_sizes.columns = ["identifier", "effect_size", "pip"]
    # SNP identifier to rsid mapping
    mapping_dct = {}
    chromosomes = [x for x in range(1, 23)]
    for chromosome in tqdm(chromosomes):
        table = pd.read_table(f"{rsid_mapping_basepath}/{chromosome}.rsid", header=None)
        for identifier, rsid in zip(table[0], table[1]):
            mapping_dct[identifier] = rsid
    # use mapping_dct to convert identifiers to rsids
    pips_and_effect_sizes["rsid"] = pips_and_effect_sizes["identifier"].map(mapping_dct)
    # save this table
    pips_and_effect_sizes.to_csv(osp.join(temp_data_dir, "all_sparsepro_effect_sizes.tsv"), sep="\t", index=False)
else:
    pips_and_effect_sizes = pd.read_csv(osp.join(temp_data_dir, "all_sparsepro_effect_sizes.tsv"), sep="\t")

# preprocessing
pips_and_effect_sizes.rename(columns={"rsid": "Name"}, inplace=True)
pips_and_effect_sizes[["chromosome", "position", "A1", "A2"]] = \
    pips_and_effect_sizes["identifier"].str.split(".", expand=True)
pips_and_effect_sizes.drop(columns=["identifier", "chromosome", "position"], inplace=True)
pips_and_effect_sizes.set_index("Name", inplace=True)

# pdb.set_trace()

# # sort pips according to decreasing pip order 
# pips_and_effect_sizes.sort_values(by="pip", ascending=False, inplace=True)
# # calculate total pip
# total_pip = pips_and_effect_sizes["pip"].sum()
# # calculate cumulative pip
# pips_and_effect_sizes["cumulative_pip"] = pips_and_effect_sizes["pip"].cumsum()

# # 95% credible set
# cs_percent = 0.95
# cs_pip = total_pip * cs_percent
# cs_95percent_table = pips_and_effect_sizes[pips_and_effect_sizes["cumulative_pip"] <= cs_pip]

# # 99% credible set 
# cs_percent = 0.99
# cs_pip = total_pip * cs_percent
# cs_99percent_table = pips_and_effect_sizes[pips_and_effect_sizes["cumulative_pip"] <= cs_pip]

# # pip > 0.01
# pips_and_effect_sizes[pips_and_effect_sizes["pip"] > 0.001].shape

for cell_type in data:
    # prepare data table 
    snp_table = data[cell_type]["original_snp_tables"]
    rsids = snp_table["Name"].tolist()
    snp_table_ = snp_table.set_index("Name")
    # prepare sparsepro data 
    sparsepro_results = pips_and_effect_sizes[pips_and_effect_sizes.index.isin(rsids)]
    sparsepro_results.columns = [f"sparsepro__{x}" for x in sparsepro_results.columns]
    # join sparsepro_results with snp table 
    joined_table = snp_table_.join(sparsepro_results, how="left")
    joined_table.reset_index(inplace=True)

    # joined_table[["A1", "A2", "sparsepro_A1", "sparsepro_A2"]]

    # sanity check (5): check if sparsepro alleles are included in A1 or A2 of alleles used for the pipeline
    sparsepro_A1 = joined_table["sparsepro__A1"].tolist()
    sparsepro_A2 = joined_table["sparsepro__A2"].tolist()
    pipeline_A1_A2 = joined_table["A1"] + "__" + joined_table["A2"]
    allele_inclusion = []
    sparsepro_alleles_exist = []
    for s_A1, s_A2, pipeline_alleles in zip(sparsepro_A1, sparsepro_A2, pipeline_A1_A2):
        try:
            if s_A1 in pipeline_alleles and s_A2 in pipeline_alleles:
                allele_inclusion.append(True)
            else:
                allele_inclusion.append(False)
            sparsepro_alleles_exist.append(True)
        except: 
            # sparsepro alleles may not exist for all SNPs
            sparsepro_alleles_exist.append(False)
            allele_inclusion.append(False)
    # account for missing sparsepro alleles
    allele_inclusion = [x for x,y in zip(allele_inclusion, sparsepro_alleles_exist) if y]
    assert all(allele_inclusion), "Sanity check (5) failed"
    joined_table["sparsepro__alleles_exist"] = sparsepro_alleles_exist



    # find the alleles that are flipped 
    pipeline_A1 = joined_table["A1"].tolist()
    pipeline_A2 = joined_table["A2"].tolist()
    is_flipped = []
    for p_A1, p_A2, s_A1, s_A2 in zip(pipeline_A1, pipeline_A2, sparsepro_A1, sparsepro_A2):
        if pd.isnull(s_A1) and pd.isnull(s_A2):
            is_flipped.append(None)
        elif p_A1 == s_A1 and p_A2 == s_A2:
            is_flipped.append(False)
        elif p_A1 == s_A2 and p_A2 == s_A1:
            is_flipped.append(True)
        else:
            pdb.set_trace()
    joined_table["sparsepro__are_pipeline_alleles_flipped"] = is_flipped

    # sanity check (6): check if snp order of the tables before and after the join are same
    assert snp_table["Name"].tolist() == joined_table["Name"].tolist(), "Sanity check (6) failed"

    # update the data table 
    data[cell_type]["original_snp_tables"] = joined_table

# save the final form of merged data structure 
with open(f"{temp_data_dir}/merged_data.v2.pkl", "wb") as f:
    pickle.dump(data, f)


print("Script finished")
