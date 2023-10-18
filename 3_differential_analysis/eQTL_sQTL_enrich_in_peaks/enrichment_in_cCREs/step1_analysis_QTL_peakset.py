import os 
import pdb 
import sys  
import warnings
import argparse
import pickle
warnings.filterwarnings("ignore")

import os.path as osp
import numpy as np
import pandas as pd
import pyranges as pr

from tqdm import tqdm
from scipy.stats import binomtest

# parse arguments from step14.jobscript.py
parser = argparse.ArgumentParser()
parser.add_argument("--peakset-name", type=str, required=True)
parser.add_argument("--qtl-name", type=str, required=True)
parser.add_argument("--out-dir", type=str, required=True)
args = parser.parse_args()

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

aggregate_fn = lambda x: ",".join([str(y) for y in x])

'''
Load hg38 genome
'''
chroms_of_interest = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

filepath = "/home/dcakma3/scratch/mdd-genome2ca/data/genomes/grch38_p13/hg38.p13.fa.fai"
hg38_p13_fai = pd.read_csv(filepath, sep="\t", header=None)

hg38_p13_fai.columns = ["chrom", "size", "offset", "line_length", "line_width"]
hg38_p13_fai = hg38_p13_fai.drop(["offset", "line_length", "line_width"], axis=1)
hg38_p13_fai = hg38_p13_fai[hg38_p13_fai["chrom"].isin(chroms_of_interest)]
hg38_p13_fai["size"] = hg38_p13_fai["size"].astype(int)
hg38_p13_size = hg38_p13_fai["size"].sum()

'''
Load peaks
'''

if args.peakset_name == "ct_M_0.05":

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

elif args.peakset_name.startswith("ct_M_0.05__splitted__"):

    ###### cell type marker peaks 
    filepath =  osp.join(prepared_peaks_dir, "broad_marker.tsv")
    table = pd.read_csv(filepath, sep="\t")
    # peaks are 1001bp, convert to 501bp
    table["Start"] = table["Start"] + 250
    table["End"] = table["End"] - 250
    # get peakName from coords
    table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
    # retain cell type of interest
    context = args.peakset_name.split("__")[-1]
    table = table[table["cell_type"] == context]
    peaks_df = table
    peaks_pr = pr.PyRanges(table, int64=True)


elif args.peakset_name == "cluster_M_0.05":

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

elif args.peakset_name.startswith("cluster_M_0.05__splitted__"):

    ###### cluster marker peaks (0.05, splitted) 
    filepath =  osp.join(prepared_peaks_dir, "subcluster_marker.tsv")
    table = pd.read_csv(filepath, sep="\t")
    # peaks are 1001bp, convert to 501bp
    table["Start"] = table["Start"] + 250
    table["End"] = table["End"] - 250
    # get peakName from coords
    table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
    # retain cell type of interest
    context = args.peakset_name.split("__")[-1]
    table = table[table["cell_type"] == context]
    peaks_df = table
    peaks_pr = pr.PyRanges(table, int64=True)


elif args.peakset_name == "ct_D_0.2":

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

elif args.peakset_name.startswith("ct_D_0.2__splitted__"):

    ###### cluster DARs (0.2) 
    filepath =  osp.join(prepared_peaks_dir, "broad_diff.tsv")
    table = pd.read_csv(filepath, sep="\t")
    # peaks are 1001bp, convert to 501bp
    table["Start"] = table["Start"] + 250
    table["End"] = table["End"] - 250
    # get peakName from coords
    table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
    # retain cell type of interest
    context = args.peakset_name.split("__")[-1]
    table = table[table["cell_type"] == context]
    peaks_df = table
    peaks_pr = pr.PyRanges(table, int64=True)

elif args.peakset_name == "cluster_D_0.2":

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

elif args.peakset_name.startswith("cluster_D_0.2__splitted__"):

    ###### cluster DARs (0.2) 
    filepath =  osp.join(prepared_peaks_dir, "subcluster_diff.tsv")
    table = pd.read_csv(filepath, sep="\t")
    # peaks are 1001bp, convert to 501bp
    table["Start"] = table["Start"] + 250
    table["End"] = table["End"] - 250
    # get peakName from coords
    table["peakName"] = table["Chromosome"] + "-" + table["Start"].astype(str) + "-" + table["End"].astype(str)
    # retain cell type of interest
    context = args.peakset_name.split("__")[-1]
    table = table[table["cell_type"] == context]
    peaks_df = table
    peaks_pr = pr.PyRanges(table, int64=True)

elif args.peakset_name == "cluster_D_0.05":

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

elif args.peakset_name.startswith("cluster_D_0.05__splitted__"):

    ###### cluster DARs (0.05) 
    filepath = osp.join(raw_data_dir, "latest_DARs/atac_fdr5.csv")
    table = pd.read_csv(filepath, sep=",")
    table = table[["gene", "cluster_id", "p_adj.loc_treatment"]]
    table[["Chromosome", "Start", "End"]] = \
        table["gene"].str.split("-", expand=True)
    table["Start"] = table["Start"].astype(int)
    table["End"] = table["End"].astype(int)
    table.rename(columns={"gene": "peakName", "cluster_id": "cluster"}, inplace=True)
    context = args.peakset_name.split("__")[-1]
    table = table[table["cluster"] == context]
    peaks_df = table
    peaks_pr = pr.PyRanges(table, int64=True)
else:
    print("Invalid peakset name")
    sys.exit()


###### eQTL and sQTL (FDR < 0.05) 
if args.qtl_name == "eQTL":
    study_names = ["BrainSeq", "CommonMind", "GTEx", "ROSMAP"]
    modality_type = "eQTL:ge"
elif args.qtl_name == "sQTL":
    study_names = ["BrainSeq", "CommonMind", "ROSMAP"]
    modality_type = "sQTL:sp"
    
modality_name = modality_type.split(":")[0]
modality_id = modality_type.split(":")[1]

rsid_gene_pairs = []
QTL_tables = []
for study_name in tqdm(study_names, desc="eQTL studies"):

    filepath = f"{study_name}.{modality_id}.prepared.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    sumstats = pd.read_csv(filepath, sep="\t")

    sumstats = sumstats[~sumstats["rsid"].isna()]
    sumstats = sumstats[~sumstats["hgnc_symbol"].isna()]
    assert (sumstats["fdr"] < 0.05).all()

    sumstats["rsid_gene"] = sumstats["rsid"] + "_" + sumstats["hgnc_symbol"]
    sumstats["rsid_gene_geneID"] = sumstats["rsid"] + "_" + sumstats["hgnc_symbol"] + "_" + sumstats["gene_id"]
    sumstats = sumstats.drop_duplicates(subset=["rsid_gene"], keep="first")

    assert sumstats["rsid_gene"].nunique() == sumstats.shape[0], "assertion failed"

    # retain pairs that were not previously available
    sumstats = sumstats[~sumstats["rsid_gene"].isin(rsid_gene_pairs)]

    # update pairs 
    rsid_gene_pairs += sumstats["rsid_gene"].tolist()
    QTL_tables.append(sumstats)

QTL_df = pd.concat(QTL_tables, axis=0)

# retain unique rsid_gene
temp = QTL_df.drop_duplicates(subset=["rsid_gene"], keep="first")
temp = temp["rsid_gene_geneID"].str.split("_", expand=True)
temp.columns = ["rsid", "gene__QTL", "gene_id__QTL"]
QTL_gene_df = temp.groupby("rsid").agg({"gene__QTL": aggregate_fn, "gene_id__QTL": aggregate_fn})

# QTL SNP table having one SNP in each row
QTL_snp_df = QTL_df.drop_duplicates(subset=["rsid_gene"], keep="first")[["rsid", "variant"]].drop_duplicates(subset="rsid", keep="first")
QTL_snp_df[["Chromosome", "Start", "A1", "A2"]] = QTL_snp_df["variant"].str.split("_", expand=True)
QTL_snp_df = QTL_snp_df.drop(columns=["A1", "A2"])
QTL_snp_df["Start"] = QTL_snp_df["Start"].astype(int)
QTL_snp_df["End"] = QTL_snp_df["Start"] + 1
QTL_snp_df = QTL_snp_df.drop(columns=["variant"])
QTL_snp_pr = pr.PyRanges(QTL_snp_df, int64=True)

# unique QTL snps
snp_df = QTL_snp_df[["rsid", "Chromosome", "Start", "End"]].drop_duplicates(subset=["rsid"], keep="first")
snp_pr = pr.PyRanges(snp_df, int64=True)

##### Analysis 

results = {}

print("QTL type:", args.qtl_name)
results["qtl_type"] = args.qtl_name
print("Peakset:", args.peakset_name)
results["peakset"] = args.peakset_name
print("Peak count:", peaks_df.shape[0])
results["peak_count"] = peaks_df.shape[0]
print("QTL (snp, QTL gene) count:", QTL_df.shape[0])
results["QTL_snp_gene_count"] = QTL_df.shape[0]
print("QTL (snp) count:", QTL_snp_df.shape[0])
results["QTL_snp_count"] = QTL_snp_df.shape[0]

snps_in_peaks_df = snp_pr.join(peaks_pr, suffix="__peak", how="left").df
snps_in_peaks_df = snps_in_peaks_df[snps_in_peaks_df["Start__peak"] != -1]
snps_in_peaks_df = snps_in_peaks_df.drop(columns=["Start__peak", "End__peak"])
__snps_in_peaks_df = snps_in_peaks_df[["rsid", "Chromosome", "Start", "End"]].drop_duplicates(subset=["rsid"], keep="first")
snps_in_peaks_pr = pr.PyRanges(__snps_in_peaks_df, int64=True)
print("QTL (snp in peak) count:", __snps_in_peaks_df.shape[0])
results["QTL_snp_in_peak_count"] = __snps_in_peaks_df.shape[0]
print("QTL (snp in peak, QTL gene) count:", snps_in_peaks_df.shape[0])
results["QTL_snp_in_peak_rsid_gene_pair_count"] = snps_in_peaks_df.shape[0]

'''
Binomial test, odds ratio, log2 fold change compared to expected
'''

num_qtls = results["QTL_snp_count"]

peakset_coverage = peaks_pr.merge().lengths().values.sum() / hg38_p13_size
results["peakset_coverage"] = peakset_coverage

num_expected_qtls_in_peak = np.ceil(peakset_coverage * num_qtls)
results["expected_qtl_count"] = num_expected_qtls_in_peak

num_qtls_in_peak = __snps_in_peaks_df.shape[0]
results["observed_qtl_count"] = num_qtls_in_peak

expected_success_rate = num_expected_qtls_in_peak / num_qtls
results["expected_success_rate"] = expected_success_rate

observed_success_rate = num_qtls_in_peak / num_qtls
results["observed_success_rate"] = observed_success_rate

log2_fold_change = np.log2(observed_success_rate / expected_success_rate)
results["log2_fold_change"] = log2_fold_change

numerator = num_qtls_in_peak / (num_qtls - num_qtls_in_peak)
denominator = num_expected_qtls_in_peak / (num_qtls - num_expected_qtls_in_peak)
OR = numerator / denominator
results["OR"] = OR

ln_OR = np.log(OR)
results["ln_OR"] = ln_OR

binomtest_result = binomtest(k = num_qtls_in_peak, n = num_qtls, p = peakset_coverage, alternative='greater')

pval = binomtest_result.pvalue
results["binomtest_pval"] = pval

neglog10_pval = -np.log10(pval)
results["neglog10_pval"] = neglog10_pval

binomtest_success_rate = binomtest_result.proportion_estimate
results["binomtest_success_rate"] = binomtest_success_rate


"""
Save results
"""
filepath = osp.join(args.out_dir, f"{args.qtl_name}__{args.peakset_name}.results.pkl")
with open(filepath, "wb") as f:
    pickle.dump(results, f)
