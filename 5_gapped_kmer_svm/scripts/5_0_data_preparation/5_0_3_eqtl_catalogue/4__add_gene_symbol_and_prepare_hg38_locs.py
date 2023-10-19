import pdb 
import os 

import os.path as osp
import numpy as np 
import pandas as pd 

from tqdm import tqdm


###### configuration 


# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# configure paths
data_id = "eqtl_catalogue"
raw_data_dir = osp.join(GKMSVM_RAW_DATA_DIR, data_id, sep = "/")
sumstats_dir = osp.join(raw_data_dir, "sumstats")
processed_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id, sep = "/")



###### Load FDR thresholded eQTL sumstats, enhance it with gene id to symbol mapping and prepare hg38 locations


study_names = ["BrainSeq", "CommonMind", "GTEx", "ROSMAP"]
modality_type = "ge"
for study_name in tqdm(study_names, desc="eQTL studies"):

    # load FDR thresholded sumstats 
    filepath = f"{study_name}.{modality_type}.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    sumstats = pd.read_csv(filepath, sep="\t")

    # set gene id as index 
    sumstats = sumstats.set_index("gene_id")

    # load gene id to symbol mapping
    filepath = f"{study_name}.{modality_type}.gene_id_mapping.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    gene_id_symbol_mapping = pd.read_csv(filepath, sep="\t")

    # remove ensembl_transcript_id and entrezgene_accession cols
    gene_id_symbol_mapping = gene_id_symbol_mapping.drop(columns=["ensembl_transcript_id", "entrezgene_accession"])

    # rename ensembl gene id to gene_id
    gene_id_symbol_mapping = gene_id_symbol_mapping.rename(columns={"ensembl_gene_id": "gene_id"})

    # set gene id as index
    gene_id_symbol_mapping = gene_id_symbol_mapping.set_index("gene_id")

    # merge sumstats with gene id to symbol mapping
    fdr_thresh_sumstats = sumstats.join(gene_id_symbol_mapping, how="left")

    # reset index 
    fdr_thresh_sumstats = fdr_thresh_sumstats.reset_index()

    # split variant column by _ to Chromosome Start A1 and A2
    fdr_thresh_sumstats[["Chromosome", "Start", "A1", "A2"]] = fdr_thresh_sumstats["variant"].str.split("_", expand=True)

    # set End column as Start + 1
    fdr_thresh_sumstats["End"] = fdr_thresh_sumstats["Start"].astype(int) + 1

    # save as tsv 
    filepath = f"{study_name}.{modality_type}.prepared.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    fdr_thresh_sumstats.to_csv(filepath, sep="\t", index=False)


###### Load FDR thresholded eQTL sumstats, enhance it with gene id to symbol mapping and prepare hg38 locations

study_names = ["BrainSeq", "CommonMind", "GTEx", "ROSMAP"]
modality_type = "sp"
for study_name in tqdm(study_names, desc="sQTL studies"):

    # load FDR thresholded sumstats 
    filepath = f"{study_name}.{modality_type}.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    sumstats = pd.read_csv(filepath, sep="\t")

    # set gene id as index 
    sumstats = sumstats.set_index("gene_id")

    # load gene id to symbol mapping
    filepath = f"{study_name}.{modality_type}.gene_id_mapping.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    gene_id_symbol_mapping = pd.read_csv(filepath, sep="\t")

    # remove ensembl_transcript_id and entrezgene_accession cols
    gene_id_symbol_mapping = gene_id_symbol_mapping.drop(columns=["ensembl_transcript_id", "entrezgene_accession"])

    # rename ensembl gene id to gene_id
    gene_id_symbol_mapping = gene_id_symbol_mapping.rename(columns={"ensembl_gene_id": "gene_id"})

    # set gene id as index
    gene_id_symbol_mapping = gene_id_symbol_mapping.set_index("gene_id")

    # merge sumstats with gene id to symbol mapping
    fdr_thresh_sumstats = sumstats.join(gene_id_symbol_mapping, how="left")

    # reset index 
    fdr_thresh_sumstats = fdr_thresh_sumstats.reset_index()

    # split variant column by _ to Chromosome Start A1 and A2
    fdr_thresh_sumstats[["Chromosome", "Start", "A1", "A2"]] = fdr_thresh_sumstats["variant"].str.split("_", expand=True)

    # set End column as Start + 1
    fdr_thresh_sumstats["End"] = fdr_thresh_sumstats["Start"].astype(int) + 1

    # save as tsv 
    filepath = f"{study_name}.{modality_type}.prepared.tsv"
    filepath = osp.join(processed_data_dir, filepath)
    fdr_thresh_sumstats.to_csv(filepath, sep="\t", index=False)

print("Script finished")