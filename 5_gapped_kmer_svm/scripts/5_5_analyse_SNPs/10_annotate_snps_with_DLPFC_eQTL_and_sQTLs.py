import pdb 
import os 
import sys 
import requests
import json
import pickle

import os.path as osp
import numpy as np 
import pandas as pd 

from tqdm import tqdm


###### configuration 

raw_data_dir = "./data/raw"

processed_data_dir = "./data/processed"
os.makedirs(processed_data_dir, exist_ok=True)

results_dir = "./results"
os.makedirs(results_dir, exist_ok=True)

# for sSNP and cdSNP tables
cols_to_retain = ['Index', 'pipeline__cell_type', 'pipeline__ssSNP', 'pipeline__tf_motif_distruption_confidence', 'snp__Name', 'snp__Chromosome_hg38', 'snp__Start_hg38', 'snp__End_hg38', 'qtilizer__gene', 'qtilizer__ensgid', 'qtilizer__colocalization', 'qtilizer__distance', 'qtilizer__tissue', 'qtilizer__p', 'qtilizer__sign_info', 'qtilizer__beta', 'qtilizer__ea', 'qtilizer__nea', 'qtilizer__source', 'qtilizer__pmid', 'qtilizer__is_sign', 'qtilizer__is_best', 'qtilizer__n_qtls', 'qtilizer__n_best', 'CISBP__match_name', 'CISBP__match_altname', 'CISBP__match_qval']

###### load ssnp related data
snp_type = "sSNP"

#### load ssnp data 
filepath = osp.join(raw_data_dir, f"{snp_type}_table.tsv")
ssnp_table = pd.read_csv(filepath, sep="\t")
ssnp_table = ssnp_table[cols_to_retain]
ssnp_table.set_index("snp__Name", inplace=True)

# create a copy for later use
ssnp_table_copy = ssnp_table.copy()

#### load eqtl catalog data and merge with data
sqtl_data_combined = []
eqtl_data_combined = []

## BrainSeq - ge
study_name = "BrainSeq"
modality_type = "ge" 

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## BrainSeq - sp
study_name = "BrainSeq"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## CommonMind - ge
study_name = "CommonMind"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## CommonMind - sp
study_name = "CommonMind"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## GTEx - ge
study_name = "GTEx"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## GTEx - sp
study_name = "GTEx"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## ROSMAP - ge
study_name = "ROSMAP"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")




## ROSMAP - sp
study_name = "ROSMAP"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
ssnp_table = ssnp_table.join(qtl_data, how="left")

# sort ssnp table by pipeline__cell_type
ssnp_table.sort_values(by="pipeline__cell_type", inplace=True)

# reset index
ssnp_table.reset_index(inplace=True)

# rename index to snp__Name
ssnp_table.rename(columns={"index": "snp__Name"}, inplace=True)

# save table to results dir 
filepath = osp.join(results_dir, f"{snp_type}_table.eqtl_catalogue.split.tsv")
ssnp_table.to_csv(filepath, sep="\t", index=True)

#### combined analysis 

## eQTL
eqtl_data_combined = pd.concat(eqtl_data_combined, axis=0)
eqtl_data_combined.sort_values(by="fdr", inplace=True)

# groupby rsid and join rows via comma
eqtl_data_combined = eqtl_data_combined.groupby("rsid").agg(aggregate)

# add eQTL__ prefix to col names
eqtl_data_combined.columns = ["eQTL__" + x for x in eqtl_data_combined.columns]

# sanity check
assert eqtl_data_combined.shape[0] == eqtl_data_combined.index.nunique()

# join with ssnp table
ssnp_table = ssnp_table_copy.join(eqtl_data_combined, how="left")

## sQTL
sqtl_data_combined = pd.concat(sqtl_data_combined, axis=0)
sqtl_data_combined.sort_values(by="fdr", inplace=True)

# groupby rsid and join rows via comma
sqtl_data_combined = sqtl_data_combined.groupby("rsid").agg(aggregate)

# add sQTL__ prefix to col names
sqtl_data_combined.columns = ["sQTL__" + x for x in sqtl_data_combined.columns]

# sanity check
assert sqtl_data_combined.shape[0] == sqtl_data_combined.index.nunique()

# join with ssnp table
ssnp_table = ssnp_table.join(sqtl_data_combined, how="left")

# sort by cell type 
ssnp_table.sort_values(by="pipeline__cell_type", inplace=True)

# reset index
ssnp_table.reset_index(inplace=True)

# rename index to snp__Name
ssnp_table.rename(columns={"index": "snp__Name"}, inplace=True)

# save to results dir 
filepath = osp.join(results_dir, f"{snp_type}_table.eqtl_catalogue.combined.tsv")
ssnp_table.to_csv(filepath, sep="\t", index=True)


###### load ssnp related data
snp_type = "cdSNP"

#### load cdsnp data 
filepath = osp.join(raw_data_dir, f"{snp_type}_table.tsv")
cdsnp_table = pd.read_csv(filepath, sep="\t")
cdsnp_table = cdsnp_table[cols_to_retain]
cdsnp_table.set_index("snp__Name", inplace=True)

# create a copy for later use
cdsnp_table_copy = cdsnp_table.copy()

#### load eqtl catalog data and merge with data
sqtl_data_combined = []
eqtl_data_combined = []

## BrainSeq - ge
study_name = "BrainSeq"
modality_type = "ge" 

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## BrainSeq - sp
study_name = "BrainSeq"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## CommonMind - ge
study_name = "CommonMind"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## CommonMind - sp
study_name = "CommonMind"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## GTEx - ge
study_name = "GTEx"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## GTEx - sp
study_name = "GTEx"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## ROSMAP - ge
study_name = "ROSMAP"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")




## ROSMAP - sp
study_name = "ROSMAP"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
cdsnp_table = cdsnp_table.join(qtl_data, how="left")

# sort ssnp table by pipeline__cell_type
cdsnp_table.sort_values(by="pipeline__cell_type", inplace=True)

# reset index
cdsnp_table.reset_index(inplace=True)

# rename index to snp__Name
cdsnp_table.rename(columns={"index": "snp__Name"}, inplace=True)

# save table to results dir 
filepath = osp.join(results_dir, f"{snp_type}_table.eqtl_catalogue.split.tsv")
cdsnp_table.to_csv(filepath, sep="\t", index=True)

#### combined analysis 

## eQTL
eqtl_data_combined = pd.concat(eqtl_data_combined, axis=0)
eqtl_data_combined.sort_values(by="fdr", inplace=True)

# groupby rsid and join rows via comma
eqtl_data_combined = eqtl_data_combined.groupby("rsid").agg(aggregate)

# add eQTL__ prefix to col names
eqtl_data_combined.columns = ["eQTL__" + x for x in eqtl_data_combined.columns]

# sanity check
assert eqtl_data_combined.shape[0] == eqtl_data_combined.index.nunique()

# join with ssnp table
cdsnp_table = cdsnp_table_copy.join(eqtl_data_combined, how="left")

## sQTL
sqtl_data_combined = pd.concat(sqtl_data_combined, axis=0)
sqtl_data_combined.sort_values(by="fdr", inplace=True)

# groupby rsid and join rows via comma
sqtl_data_combined = sqtl_data_combined.groupby("rsid").agg(aggregate)

# add sQTL__ prefix to col names
sqtl_data_combined.columns = ["sQTL__" + x for x in sqtl_data_combined.columns]

# sanity check
assert sqtl_data_combined.shape[0] == sqtl_data_combined.index.nunique()

# join with ssnp table
cdsnp_table = cdsnp_table.join(sqtl_data_combined, how="left")

# sort by cell type 
cdsnp_table.sort_values(by="pipeline__cell_type", inplace=True)

# reset index
cdsnp_table.reset_index(inplace=True)

# rename index to snp__Name
cdsnp_table.rename(columns={"index": "snp__Name"}, inplace=True)

# save to results dir 
filepath = osp.join(results_dir, f"{snp_type}_table.eqtl_catalogue.combined.tsv")
cdsnp_table.to_csv(filepath, sep="\t", index=True)












###### load ssnp related data
snp_type = "MDD"

#### load ssnp data 
filepath = osp.join(raw_data_dir, f"{snp_type}_relevant_snps.tsv")
mdd_snp_table = pd.read_csv(filepath, sep="\t")
mdd_snp_table.set_index("Name", inplace=True)

# create a copy for later use
mdd_snp_table_copy = mdd_snp_table.copy()

#### load eqtl catalog data and merge with data
sqtl_data_combined = []
eqtl_data_combined = []

## BrainSeq - ge
study_name = "BrainSeq"
modality_type = "ge" 

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## BrainSeq - sp
study_name = "BrainSeq"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## CommonMind - ge
study_name = "CommonMind"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## CommonMind - sp
study_name = "CommonMind"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## GTEx - ge
study_name = "GTEx"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## GTEx - sp
study_name = "GTEx"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## ROSMAP - ge
study_name = "ROSMAP"
modality_type = "ge"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")




## ROSMAP - sp
study_name = "ROSMAP"
modality_type = "sp"

# load results 
filepath = f"{snp_type}.{study_name}.{modality_type}.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data = pd.read_csv(filepath, sep="\t")
qtl_data.set_index("gene_id", inplace=True)
qtl_data["study_name"] = study_name
qtl_data["modality_type"] = modality_type

# load gene symbol mapping 
filepath = f"{snp_type}.{study_name}.{modality_type}.gene_symbol_mapping.tsv"
filepath = osp.join(processed_data_dir, filepath)
qtl_data_mapping = pd.read_csv(filepath, sep="\t")

# prepare gene symbol mapping table
qtl_data_mapping.drop("ensembl_transcript_id", axis=1, inplace=True) 
qtl_data_mapping.drop_duplicates(subset="ensembl_gene_id", keep="first", inplace=True)
qtl_data_mapping.set_index("ensembl_gene_id", inplace=True)

# add gene symbol to qtl_data
qtl_data = qtl_data.join(qtl_data_mapping, how="left")
qtl_data.reset_index(inplace=True)
qtl_data.rename(columns={"index": "gene_id"}, inplace=True)

if modality_type == "ge":
    eqtl_data_combined.append(qtl_data)
elif modality_type == "sp":
    sqtl_data_combined.append(qtl_data)
else:
    raise ValueError(f"modality_type {modality_type} not recognized for study {study_name}")

# sort by increasing fdr
qtl_data.sort_values(by="fdr", inplace=True)

# group by rsid and join rows via comma
def aggregate(x):
    return ",".join([str(y) if not pd.isna(y) else "" for y in x])
qtl_data = qtl_data.groupby("rsid").agg(aggregate)

# sanity check 
assert qtl_data.shape[0] == qtl_data.index.nunique()

# add prefix to col names
qtl_data.columns = [f"{study_name}__{modality_type}__" + x for x in qtl_data.columns]

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(qtl_data, how="left")

# reset index
mdd_snp_table.reset_index(inplace=True)

# rename index to snp__Name
mdd_snp_table.rename(columns={"index": "Name"}, inplace=True)

# save table to results dir 
filepath = osp.join(results_dir, f"{snp_type}_snp_table.eqtl_catalogue.split.tsv")
mdd_snp_table.to_csv(filepath, sep="\t", index=True)

#### combined analysis 

## eQTL
eqtl_data_combined = pd.concat(eqtl_data_combined, axis=0)
eqtl_data_combined.sort_values(by="fdr", inplace=True)

# groupby rsid and join rows via comma
eqtl_data_combined = eqtl_data_combined.groupby("rsid").agg(aggregate)

# add eQTL__ prefix to col names
eqtl_data_combined.columns = ["eQTL__" + x for x in eqtl_data_combined.columns]

# sanity check
assert eqtl_data_combined.shape[0] == eqtl_data_combined.index.nunique()

# join with ssnp table
mdd_snp_table = mdd_snp_table_copy.join(eqtl_data_combined, how="left")

## sQTL
sqtl_data_combined = pd.concat(sqtl_data_combined, axis=0)
sqtl_data_combined.sort_values(by="fdr", inplace=True)

# groupby rsid and join rows via comma
sqtl_data_combined = sqtl_data_combined.groupby("rsid").agg(aggregate)

# add sQTL__ prefix to col names
sqtl_data_combined.columns = ["sQTL__" + x for x in sqtl_data_combined.columns]

# sanity check
assert sqtl_data_combined.shape[0] == sqtl_data_combined.index.nunique()

# join with ssnp table
mdd_snp_table = mdd_snp_table.join(sqtl_data_combined, how="left")

# reset index
mdd_snp_table.reset_index(inplace=True)

# rename index to snp__Name
mdd_snp_table.rename(columns={"index": "Name"}, inplace=True)

# save to results dir 
filepath = osp.join(results_dir, f"{snp_type}_snp_table.eqtl_catalogue.combined.tsv")
mdd_snp_table.to_csv(filepath, sep="\t", index=True)

print("Script finished")
