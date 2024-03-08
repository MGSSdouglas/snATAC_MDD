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

# cell types to analyze
cell_types = ["Ast1", "Ast2", "Ast3", "Ast4", "End2", "ExN1", \
    "ExN1_L23", "ExN1_L24", "ExN1_L46", "ExN1_L56", "ExN2", "ExN2_L23", \
    "ExN2_L46", "ExN2_L56", "ExN3_L46", "ExN3_L56", "ExN4_L56", "In_LAMP5", \
    "InN3", "In_PV", "In_SST", "In_VIP", "Mic1", "Mic2", "End1", "Oli1", "Oli2", \
    "Oli3", "Oli4", "Oli5", "Oli6", "Oli7", "OPC1", "OPC2", "OPC3", "OPC4"]

broad2sub = {
    "Ast": ["Ast1", "Ast2", "Ast3", "Ast4"],
    "End": ["End2"],
    "ExN": ["ExN1", "ExN1_L23", "ExN1_L24", "ExN1_L46", "ExN1_L56", "ExN2", "ExN2_L23", "ExN2_L46", "ExN2_L56", "ExN3_L46", "ExN3_L56", "ExN4_L56"],
    "InN": ["In_LAMP5", "InN3", "In_PV", "In_SST", "In_VIP"],
    "Mic": ["Mic1", "Mic2", "End1"],
    "Oli": ["Oli1", "Oli2", "Oli3", "Oli4", "Oli5", "Oli6", "Oli7"],
    "OPC": ["OPC1", "OPC2", "OPC3", "OPC4"]
}

broad_clusters = []
for x in cell_types:
    if x in broad2sub["Ast"]:
        broad_clusters.append("Ast")
    elif x in broad2sub["End"]:
        broad_clusters.append("End")
    elif x in broad2sub["ExN"]:
        broad_clusters.append("ExN")
    elif x in broad2sub["InN"]:
        broad_clusters.append("InN")
    elif x in broad2sub["Mic"]:
        broad_clusters.append("Mic")
    elif x in broad2sub["Oli"]:
        broad_clusters.append("Oli")
    elif x in broad2sub["OPC"]:
        broad_clusters.append("OPC")
    else:
        raise ValueError("Cell type {} not found in broad2sub dictionary".format(x))

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


# load merged data structures 
with open(f"{temp_data_dir}/merged_data.v2.pkl", "rb") as f:
    data = pickle.load(f)

# Steps to be performed per cell type
# 1. Calculate null importance score distributions using snps and null sequences
#   fit t-distiribution and find 95% confidence interval bounds according to the fitted distribution
#   a. ISM
#   b. gkmexplain
#   c. deltaSVM
# 2. Find statistically significant SNPs according to each score type
#   a. ISM
#   b. gkmexplain
#   c. deltaSVM
# 3. Overlap significant SNPs across score types and find those that 
#   are statistically significant across all three score types
def find_ssSNP_status(scores, lower_ci, upper_ci):
    return [True if score > upper_ci or score < lower_ci else False for score in scores]

for cell_type in data:
    plot_base_dir = f"{output_plot_dir}/{cell_type}"
    os.makedirs(plot_base_dir, exist_ok=True)

    """
    Calculate null score distributions
    """
    score_types = ["ISM", "gkmexplain", "deltaSVM"]
    for score_type in score_types:

        data[cell_type][f"null_{score_type}_ve_score_related"] = {}
        # get null scores
        null_scores = data[cell_type][f"shuffled_{score_type}_ve_scores"]
        # plot histogram of null scores 
        plt.figure()
        plt.hist(null_scores, bins=20, color="dimgray")
        plt.xlabel(f"Null {score_type} A1 - A2 scores")
        plt.ylabel("# of observations")
        plt.savefig(f"./{plot_base_dir}/null_{score_type}_ve_score_distribution.png")
        plt.close()
        #fit t distribution to the data 
        t_df, t_loc, t_scale = stats.t.fit(null_scores)
        t = stats.t(df=t_df, loc=t_loc, scale=t_scale)
        data[cell_type][f"null_{score_type}_ve_score_related"]["fitted_t_df"] = t_df
        data[cell_type][f"null_{score_type}_ve_score_related"]["fitted_t_loc"] = t_loc
        data[cell_type][f"null_{score_type}_ve_score_related"]["fitted_t_scale"] = t_scale
        data[cell_type][f"null_{score_type}_ve_score_related"]["fitted_t"] = t
        # apply kstest to null VE scores with t distribution and normal distribution
        t_dist_ks_test_result = stats.kstest(null_scores, "t", \
            args=(t_df, t_loc, t_scale), alternative="two_sided", \
            N=null_scores.shape[0])
        data[cell_type][f"null_{score_type}_ve_score_related"]["fitted_t_KS_p_val"] = t_dist_ks_test_result.pvalue
        data[cell_type][f"null_{score_type}_ve_score_related"]["fitted_t_KS_statistic"] = t_dist_ks_test_result.statistic
        print(f"[INFO, {cell_type}] null {score_type} ve score KS test: p-val, statistic = {t_dist_ks_test_result.pvalue}, {t_dist_ks_test_result.statistic}")
        # find 95% confidence interval bounds
        lower_ci, upper_ci = t.interval(0.95)
        data[cell_type][f"null_{score_type}_ve_score_related"]["95%_CI_upper"] = upper_ci
        data[cell_type][f"null_{score_type}_ve_score_related"]["95%_CI_lower"] = lower_ci
        # plot null score distribution, overlay t distribution and CI
        x = np.arange(null_scores.min(), null_scores.max(), 0.001)
        plt.figure()
        plt.hist(null_scores, bins=1000, density=True, color="dimgray")
        plt.plot(x, t.pdf(x), "r-", linewidth=2)
        plt.axvline(lower_ci, color="y", linestyle="-")
        plt.axvline(upper_ci, color="y", linestyle="-")
        plt.xlabel(f"Null {score_type} A1 - A2 scores")
        plt.ylabel("Density")
        plt.title(f"t-distr, two-sided KS test p-val = {t_dist_ks_test_result.pvalue}, statistic = {t_dist_ks_test_result.statistic}")
        plt.savefig(f"./{plot_base_dir}/null_{score_type}_ve_score_distribution.overlay_t_and_CI.png")
        plt.close()
        # find statistically significant SNPs according to the current score
        snp_table = data[cell_type]["original_snp_tables"]
        snp_table[f"{score_type}_ssSNP"] = find_ssSNP_status(snp_table[f"{score_type}_ve_score"], lower_ci, upper_ci)
        print(f"[INFO, {cell_type}] the number of ssSNPs according to {score_type}: {snp_table[f'{score_type}_ssSNP'].sum()}")
    
    # find statistically significant SNPs across all three scores
    snp_table["ssSNP"] = snp_table[f"gkmexplain_ssSNP"] & snp_table[f"ISM_ssSNP"] & snp_table[f"deltaSVM_ssSNP"]
    print(F"[INFO, {cell_type}] the number of ssSNPs: {snp_table['ssSNP'].sum()}")

    # make changes to the original table 
    data[cell_type]["original_snp_tables"] = snp_table

    # add cell type to the snp table 
    data[cell_type]["original_snp_tables"]["analysis_cell_type"] = cell_type

# save modified data structure 
with open(f"{temp_data_dir}/merged_data.v3.pkl", "wb") as f:
    pickle.dump(data, f)

# flatten snp_tables across cell types for eQTL analysis
tables = []
for cell_type in data:
    tables.append(data[cell_type]["original_snp_tables"])
snp_table = pd.concat(tables, axis=0)

# save snp table
snp_table.to_csv(f"{temp_data_dir}/merged_data.v3.ct_independent.csv", index=False)


print("Script finished")
