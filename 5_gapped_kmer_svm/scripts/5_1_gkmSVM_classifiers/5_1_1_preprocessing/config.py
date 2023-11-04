import os
import pdb
import pickle

import os.path as osp
import pandas as pd
import numpy as np

from Bio.SeqUtils import GC


"""
This script contains the configuration and data preparation script for the
compute job for generating cluster-specific negative examples 
not overlapping with OCRs of the cluster
"""

## Configuration

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_ENV_DIR = os.environ.get("GKMSVM_ENV_DIR")
GKMSVM_SCRIPTS_DIR = os.environ.get("GKMSVM_SCRIPTS_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# create temporary path for current analysis
temp_storage_basepath = osp.join(GKMSVM_TMP_DIR, "5_1_gkmSVM_classifiers")
os.makedirs(temp_storage_basepath, exist_ok=True)

# load prepared hg38 genome
genome_path = osp.join(GKMSVM_RAW_DATA_DIR, "hg38/genomes/grch38_p13.dict.pkl")
with open(genome_path, "rb") as f:
    genome = pickle.load(f)

# generate candidate 1001-bp negative examples by 
# tiling hg38 genome 1001bp at a time with 50bp stride
seqlen = 1001
stride = 50
compute_GC = lambda sequence: GC(sequence) / 100
genome_df = dict(Chromosome=[], Start=[], End=[], GC_1001bp=[])
for chrom, chrom_seq in genome.items():
    for i in range(0, len(chrom_seq)-seqlen, stride):
        seq = chrom_seq[i:i+seqlen]
        if "N" in seq or "n" in seq:
            continue
        else:
            genome_df["Start"].append(i)
            genome_df["End"].append(i+seqlen-1)
            genome_df["Chromosome"].append(chrom)
            genome_df["GC_1001bp"].append(compute_GC(seq))
genome_df = pd.DataFrame(genome_df)
save_path = osp.join(GKMSVM_TMP_DIR, "5_1_gkmSVM_classifiers/candidate_negative_examples.tsv")
genome_df.to_csv(save_path, sep="\t")

## prepare jobs

# source filepaths
ocr_dir = osp.join(GKMSVM_RAW_DATA_DIR, "peaks/cluster__OCRs")
clusters  = [x.split("_OCR.tsv") for x in os.listdir(ocr_dir)]
candidate_negative_path = osp.join(GKMSVM_TMP_DIR, "5_1_gkmSVM_classifiers/candidate_negative_examples.tsv")

# target filepaths
temp_storage_basepath = osp.join(GKMSVM_TMP_DIR, "5_1_gkmSVM_classifiers/5_1_1_preprocessing")
os.makedirs(temp_storage_basepath, exist_ok=True)

# const
fold_split_filepath = "./cv_splits.txt"
num_gc_bins = 20

# job params
jobs = []
for cluster in clusters:
    jobs.append({"Name": cluster,
                "cluster": cluster,
                "ocr_dir": ocr_dir,
                "candidate_neg_path": candidate_negative_path,
                "target_dir": temp_storage_basepath,
                "genome_path": genome_path,
                "fold_split_file":fold_split_filepath,
                "num_gc_bins":num_gc_bins})
