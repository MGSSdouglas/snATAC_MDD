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

"""
This script loads the FDR thresholded eQTL and sQTL sumstats and saves the unique gene ids as a tsv file.
"""

###### configuration 

raw_data_dir = "./data/raw"

processed_data_dir = "./data/processed"
os.makedirs(processed_data_dir, exist_ok=True)

results_dir = "./results"
os.makedirs(results_dir, exist_ok=True)


###### Load FDR thresholded eQTL sumstats and save as tsv 

study_names = ["BrainSeq", "CommonMind", "GTEx", "ROSMAP"]
modality_type = "ge"
unique_gene_ids = []
for study_name in tqdm(study_names, desc="eQTL studies"):

    # load FDR thresholded sumstats 
    filepath = f"{study_name}.{modality_type}.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    sumstats = pd.read_csv(filepath, sep="\t")

    # get unique genes and convert to dataframe
    gene_ids = sumstats["gene_id"].unique().tolist()
    gene_id_df = pd.DataFrame(dict(gene_id=gene_ids))

    # save as tsv 
    filepath = f"{study_name}.{modality_type}.gene_ids.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    gene_id_df.to_csv(filepath, sep="\t", index=False)

    unique_gene_ids += [x for x in gene_ids if x not in unique_gene_ids]

print(f"Number of unique genes: {len(unique_gene_ids)}")


###### Load FDR thresholded sQTL sumstats and save as tsv 

study_names = ["BrainSeq", "CommonMind", "GTEx", "ROSMAP"]
modality_type = "sp"
unique_gene_ids = []
for study_name in tqdm(study_names, desc="sQTL studies"):

    # load FDR thresholded sumstats 
    filepath = f"{study_name}.{modality_type}.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    sumstats = pd.read_csv(filepath, sep="\t")

    # get unique rsids and convert to dataframe
    gene_ids = sumstats["gene_id"].unique().tolist()
    gene_id_df = pd.DataFrame(dict(gene_id=gene_ids))

    # save as tsv 
    filepath = f"{study_name}.{modality_type}.gene_ids.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    gene_id_df.to_csv(filepath, sep="\t", index=False)

    unique_gene_ids += [x for x in gene_ids if x not in unique_gene_ids]

print(f"Number of unique genes: {len(unique_gene_ids)}")


print("Script finished")
