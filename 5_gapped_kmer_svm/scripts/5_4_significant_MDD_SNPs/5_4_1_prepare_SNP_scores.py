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


"""
Configuration
"""

# cross validation fold count
num_folds = 5
null_sequences_count = 10

seq_len = "201bp"
seq_types = ["original", "shuffled"]

effect_allele = "a1"
non_effect_allele = "a2"

# broad cluster data paths 
analysis_type = "broad"

broad_gkmexplain_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain"
broad_deltaSVM_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/deltaSVM"
broad_ISM_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/ISM"
broad_snp_info_path = f"/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/{analysis_type}"

# subcluster data paths 
analysis_type = "subcluster"

subcluster_gkmexplain_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain"
subcluster_ISM_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/ISM"
subcluster_deltaSVM_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/deltaSVM"
subcluster_snp_info_path = f"/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/{analysis_type}"

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"
os.makedirs(temp_data_dir, exist_ok=True)

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
os.makedirs(output_data_dir, exist_ok=True)


"""
Load broad cluster data
"""
print("[INFO] Loading broad cluster data")
broad_cluster_cell_types_to_omit = []

''' Original SNP tables '''
broad__original_snp_tables = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    try:
        broad__original_snp_tables[cell_type] = \
            pd.read_csv(f"{broad_snp_info_path}/{broad_cluster}/{cell_type}/snps_in_peaks.tsv", sep="\t")
    except:
        print(f"[WARNING] (broad) {cell_type} SNP table does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
        continue

''' Create Shuffled SNP tables from original SNP tables '''
broad__shuffled_snp_tables = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    try:
        broad__shuffled_snp_tables[cell_type] = \
            pd.read_csv(f"{broad_snp_info_path}/{broad_cluster}/{cell_type}/snps_in_peaks.tsv", sep="\t")
        # repeat each row null_sequences_count times without changing the order
        broad__shuffled_snp_tables[cell_type] = \
            broad__shuffled_snp_tables[cell_type].reindex(
                broad__shuffled_snp_tables[cell_type].index.repeat(null_sequences_count)
            )
    except:
        print(f"[WARNING] (broad) {cell_type} SNP table does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
        continue

''' Original 201bp SNP sequences '''
broad__original_201bp_snp_sequences_a1_allele = {}
broad__original_201bp_snp_sequences_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # 201bp a1 allele original fasta path 
    a1_allele_fasta = osp.join(*[
        broad_snp_info_path,
        broad_cluster,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a1.fa"
    ])
    if not osp.exists(a1_allele_fasta):
        print(f"[WARNING] (Original, broad) {cell_type} A1 allele sequence does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__original_201bp_snp_sequences_a1_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a1_allele_fasta, "fasta")]
    # 201bp a2 allele original fasta path
    a2_allele_fasta = osp.join(*[
        broad_snp_info_path,
        broad_cluster,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a2.fa"
    ])
    if not osp.exists(a2_allele_fasta):
        print(f"[WARNING] (Original, broad) {cell_type} A2 allele sequence does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__original_201bp_snp_sequences_a2_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a2_allele_fasta, "fasta")]

''' Shuffled 201bp SNP sequences '''
broad__shuffled_201bp_snp_sequences_a1_allele = {}
broad__shuffled_201bp_snp_sequences_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # 201bp a1 allele original fasta path 
    a1_allele_fasta = osp.join(*[
        broad_snp_info_path,
        broad_cluster,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a1.shuffled.fa"
    ])
    if not osp.exists(a1_allele_fasta):
        print(f"[WARNING] (Original, broad) {cell_type} A1 allele sequence does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__shuffled_201bp_snp_sequences_a1_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a1_allele_fasta, "fasta")]
    # 201bp a2 allele original fasta path
    a2_allele_fasta = osp.join(*[
        broad_snp_info_path,
        broad_cluster,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a2.shuffled.fa"
    ])
    if not osp.exists(a2_allele_fasta):
        print(f"[WARNING] (Original, broad) {cell_type} A2 allele sequence does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__shuffled_201bp_snp_sequences_a2_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a2_allele_fasta, "fasta")]

''' Original sequence gkmexplain variant effect scores '''
broad__original_gkmexplain_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    ve_scores = osp.join(*[
        broad_gkmexplain_filepath_base,
        broad_cluster,
        cell_type,
        f"normalized.gkmexplain_variant_effect_score.{seq_len}.original.imp_score.npy"
    ])
    if not osp.exists(ve_scores):
        print(f"[WARNING] (Original, broad) {cell_type} scores does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
        continue
    broad__original_gkmexplain_variant_effect_scores[cell_type] = np.load(ve_scores)

''' Shuffled sequence gkmexplain variant effect scores '''
broad__shuffled_gkmexplain_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    ve_scores = osp.join(*[
        broad_gkmexplain_filepath_base,
        broad_cluster,
        cell_type,
        f"normalized.gkmexplain_variant_effect_score.{seq_len}.shuffled.imp_score.npy"
    ])
    if not osp.exists(ve_scores):
        print(f"[WARNING] (Shuffled, broad) {cell_type} scores does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
        continue
    broad__shuffled_gkmexplain_variant_effect_scores[cell_type] = np.load(ve_scores)

''' Original sequence allele specific 201bp gkmexplain scores '''
broad__original_201bp_gkmexplain_scores_a1_allele = {}
broad__original_201bp_gkmexplain_scores_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # load a1 scores 
    a1_scores = osp.join(*[
        broad_gkmexplain_filepath_base,
        broad_cluster,
        cell_type,
        f"snps.hg38_information.{seq_len}.a1.original.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a1_scores):
        print(f"[WARNING] (Original, broad) {cell_type} A1 allele gkmexplain scores does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__original_201bp_gkmexplain_scores_a1_allele[cell_type] = np.load(a1_scores)
    # load a2 scores
    a2_scores = osp.join(*[
        broad_gkmexplain_filepath_base,
        broad_cluster,
        cell_type,
        f"snps.hg38_information.{seq_len}.a2.original.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a2_scores):
        print(f"[WARNING] (Original, broad) {cell_type} A2 allele gkmexplain scores does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__original_201bp_gkmexplain_scores_a2_allele[cell_type] = np.load(a2_scores)

''' Shuffled sequence allele specific 201bp gkmexplain scores '''
broad__shuffled_201bp_gkmexplain_scores_a1_allele = {}
broad__shuffled_201bp_gkmexplain_scores_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # load a1 scores 
    a1_scores = osp.join(*[
        broad_gkmexplain_filepath_base,
        broad_cluster,
        cell_type,
        f"snps.hg38_information.{seq_len}.a1.shuffled.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a1_scores):
        print(f"[WARNING] (Shuffled, broad) {cell_type} A1 allele gkmexplain scores does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__shuffled_201bp_gkmexplain_scores_a1_allele[cell_type] = np.load(a1_scores)
    # load a2 scores
    a2_scores = osp.join(*[
        broad_gkmexplain_filepath_base,
        broad_cluster,
        cell_type,
        f"snps.hg38_information.{seq_len}.a2.shuffled.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a2_scores):
        print(f"[WARNING] (Shuffled, broad) {cell_type} A2 allele gkmexplain scores does not exist")
        broad_cluster_cell_types_to_omit += [cell_type]
    else:
        broad__shuffled_201bp_gkmexplain_scores_a2_allele[cell_type] = np.load(a2_scores)

''' Original sequence deltaSVM scores '''
broad__original_deltaSVM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        broad_deltaSVM_filepath_base,
        broad_cluster,
        cell_type, 
        f"snps.hg38_information.{seq_len}.original.{effect_allele}_{non_effect_allele}.deltasvm.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Original, broad) {} deltaSVM scores do not exist".format(cell_type))
        broad_cluster_cell_types_to_omit += [cell_type]
        continue
    broad__original_deltaSVM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_deltasvm_score"].tolist())

''' Shuffled sequence deltaSVM scores '''
broad__shuffled_deltaSVM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        broad_deltaSVM_filepath_base,
        broad_cluster,
        cell_type, 
        f"snps.hg38_information.{seq_len}.shuffled.{effect_allele}_{non_effect_allele}.deltasvm.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Shuffled, broad) {} deltaSVM scores do not exist".format(cell_type))
        broad_cluster_cell_types_to_omit += [cell_type]
        continue
    broad__shuffled_deltaSVM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_deltasvm_score"].tolist())

''' Original sequence ISM scores '''
broad__original_ISM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        broad_ISM_filepath_base,
        broad_cluster,
        cell_type, 
        f"snps.hg38_information.{seq_len}.original.{effect_allele}_{non_effect_allele}.ism.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Original, broad) {} ISM scores do not exist".format(cell_type))
        broad_cluster_cell_types_to_omit += [cell_type]
        continue
    broad__original_ISM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_ism_score"].tolist())

''' Shuffled sequence ISM scores '''
broad__shuffled_ISM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        broad_ISM_filepath_base,
        broad_cluster,
        cell_type, 
        f"snps.hg38_information.{seq_len}.shuffled.{effect_allele}_{non_effect_allele}.ism.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Shuffled, broad) {} ISM scores do not exist".format(cell_type))
        broad_cluster_cell_types_to_omit += [cell_type]
        continue
    broad__shuffled_ISM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_ism_score"].tolist())

broad_cluster_cell_types_to_omit = list(set(broad_cluster_cell_types_to_omit))

broad_cluster_data = {
	"original_snp_tables": broad__original_snp_tables, # done
	"shuffled_snp_tables": broad__shuffled_snp_tables, # done
	"original_201bp_a1_allele_sequence": broad__original_201bp_snp_sequences_a1_allele, # done
	"original_201bp_a2_allele_sequence": broad__original_201bp_snp_sequences_a2_allele, # done
	"shuffled_201bp_a1_allele_sequence": broad__shuffled_201bp_snp_sequences_a1_allele, # done
	"shuffled_201bp_a2_allele_sequence": broad__shuffled_201bp_snp_sequences_a2_allele, # done
	"original_gkmexplain_ve_scores": broad__original_gkmexplain_variant_effect_scores, # done
	"shuffled_gkmexplain_ve_scores": broad__shuffled_gkmexplain_variant_effect_scores, # done
	"original_201bp_gkmexplain_a1_allele_scores": broad__original_201bp_gkmexplain_scores_a1_allele, # done
	"original_201bp_gkmexplain_a2_allele_scores": broad__original_201bp_gkmexplain_scores_a2_allele, # done
	"shuffled_201bp_gkmexplain_a1_allele_scores": broad__shuffled_201bp_gkmexplain_scores_a1_allele, # done
	"shuffled_201bp_gkmexplain_a2_allele_scores": broad__shuffled_201bp_gkmexplain_scores_a2_allele, # dones
	"original_deltaSVM_ve_scores": broad__original_deltaSVM_variant_effect_scores, # done
	"shuffled_deltaSVM_ve_scores": broad__shuffled_deltaSVM_variant_effect_scores, # done
	"original_ISM_ve_scores": broad__original_ISM_variant_effect_scores, # done
	"shuffled_ISM_ve_scores": broad__shuffled_ISM_variant_effect_scores, # done
    "cell_types_to_omit": broad_cluster_cell_types_to_omit
}

''' Sanity check '''
broad__cell_types_to_retain = list(set(cell_types) - set(broad_cluster_data['cell_types_to_omit']))
print("[INFO] The # of cell types to retain for broad cluster data is {}".format(len(broad__cell_types_to_retain)))
for cell_type in broad__cell_types_to_retain:

    # check if these cell types are available within the data structure
    original_data = []
    shuffled_data = []
    for outer_key, data in broad_cluster_data.items():
        if outer_key == "cell_types_to_omit":
            pass
        elif cell_type not in data:
            print(f"[ERROR] {cell_type} is not available in {outer_key} data")
        elif data[cell_type] is None:
            print(f"[ERROR] {cell_type} is available in {outer_key} data but empty!")
        elif "original" in outer_key:
            original_data += [data[cell_type]]
        elif "shuffled" in outer_key:
            shuffled_data += [data[cell_type]]
        elif "cell_types_to_omit" in outer_key:
            pass
        else:
            pdb.set_trace()
            print("[ERROR] Sanity check (edge case)")

    # check if dimensions of the original data are coherent
    original_data_shapes = [len(x) for x in original_data]
    if len(set(original_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR] Sanity check (dimensions of the original data are not coherent)")
    else:
        print(f"[INFO] {cell_type} original data dimension is {original_data_shapes[0]}")
    
    # check if dimensions of the shuffled data are coherent
    shuffled_data_shapes = [len(x) for x in shuffled_data]
    if len(set(shuffled_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR] Sanity check (dimensions of the shuffled data are not coherent)")
    else:
        print(f"[INFO] {cell_type} shuffled data dimension is {shuffled_data_shapes[0]}")



"""
Load subcluster data
"""
print("[INFO] Loading subcluster data")
subcluster_cell_types_to_omit = []

''' Original SNP tables '''
subcluster__original_snp_tables = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    try:
        subcluster__original_snp_tables[cell_type] = \
            pd.read_csv(f"{subcluster_snp_info_path}/{cell_type}/snps_in_peaks.tsv", sep="\t")
    except:
        print(f"[WARNING] (Original, subcluster) {cell_type} SNP table does not exist")
        subcluster_cell_types_to_omit += [cell_type]
        continue

''' Create Shuffled SNP tables from original SNP tables '''
subcluster__shuffled_snp_tables = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    try:
        subcluster__shuffled_snp_tables[cell_type] = \
            pd.read_csv(f"{subcluster_snp_info_path}/{cell_type}/snps_in_peaks.tsv", sep="\t")
        # repeat each row null_sequences_count times without changing the order
        subcluster__shuffled_snp_tables[cell_type] = \
            subcluster__shuffled_snp_tables[cell_type].reindex(
                subcluster__shuffled_snp_tables[cell_type].index.repeat(null_sequences_count)
            )
    except:
        print(f"[WARNING] (Shuffled, subcluster) {cell_type} SNP table could not be created. Please check the original snp table")
        subcluster_cell_types_to_omit += [cell_type]
        continue

''' Original 201bp SNP sequences '''
subcluster__original_201bp_snp_sequences_a1_allele = {}
subcluster__original_201bp_snp_sequences_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # 201bp a1 allele original fasta path 
    a1_allele_fasta = osp.join(*[
        subcluster_snp_info_path,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a1.fa"
    ])
    if not osp.exists(a1_allele_fasta):
        print(f"[WARNING] (Original, subcluster) {cell_type} A1 allele sequence does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__original_201bp_snp_sequences_a1_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a1_allele_fasta, "fasta")]
    # 201bp a2 allele original fasta path
    a2_allele_fasta = osp.join(*[
        subcluster_snp_info_path,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a2.fa"
    ])
    if not osp.exists(a2_allele_fasta):
        print(f"[WARNING] (Original, subcluster) {cell_type} A2 allele sequence does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__original_201bp_snp_sequences_a2_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a2_allele_fasta, "fasta")]

''' Shuffled 201bp SNP sequences '''
subcluster__shuffled_201bp_snp_sequences_a1_allele = {}
subcluster__shuffled_201bp_snp_sequences_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # 201bp a1 allele original fasta path 
    a1_allele_fasta = osp.join(*[
        subcluster_snp_info_path,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a1.shuffled.fa"
    ])
    if not osp.exists(a1_allele_fasta):
        print(f"[WARNING] (Shuffled, subcluster) {cell_type} A1 allele sequence does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__shuffled_201bp_snp_sequences_a1_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a1_allele_fasta, "fasta")]
    # 201bp a2 allele original fasta path
    a2_allele_fasta = osp.join(*[
        subcluster_snp_info_path,
        cell_type,
        f"{seq_len}/snps.hg38_information.{seq_len}.a2.shuffled.fa"
    ])
    if not osp.exists(a2_allele_fasta):
        print(f"[WARNING] (Shuffled, subcluster) {cell_type} A2 allele sequence does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__shuffled_201bp_snp_sequences_a2_allele[cell_type] = \
            [str(x.seq) for x in SeqIO.parse(a2_allele_fasta, "fasta")]

''' Original sequence gkmexplain variant effect scores '''
subcluster__original_gkmexplain_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    ve_scores = osp.join(*[
        subcluster_gkmexplain_filepath_base,
        cell_type,
        f"normalized.gkmexplain_variant_effect_score.{seq_len}.original.imp_score.npy"
    ])
    if not osp.exists(ve_scores):
        print(f"[WARNING] (Original, subcluster) {cell_type} scores does not exist")
        subcluster_cell_types_to_omit += [cell_type]
        continue
    subcluster__original_gkmexplain_variant_effect_scores[cell_type] = np.load(ve_scores)

''' Shuffled sequence gkmexplain variant effect scores '''
subcluster__shuffled_gkmexplain_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    ve_scores = osp.join(*[
        subcluster_gkmexplain_filepath_base,
        cell_type,
        f"normalized.gkmexplain_variant_effect_score.{seq_len}.shuffled.imp_score.npy"
    ])
    if not osp.exists(ve_scores):
        print(f"[WARNING] (Shuffled, subcluster) {cell_type} scores does not exist")
        subcluster_cell_types_to_omit += [cell_type]
        continue
    subcluster__shuffled_gkmexplain_variant_effect_scores[cell_type] = np.load(ve_scores)

''' Original sequence allele specific 201bp gkmexplain scores '''
subcluster__original_201bp_gkmexplain_scores_a1_allele = {}
subcluster__original_201bp_gkmexplain_scores_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # load a1 scores 
    a1_scores = osp.join(*[
        subcluster_gkmexplain_filepath_base,
        cell_type,
        f"snps.hg38_information.{seq_len}.a1.original.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a1_scores):
        print(f"[WARNING] (Original, subcluster) {cell_type} A1 allele gkmexplain scores does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__original_201bp_gkmexplain_scores_a1_allele[cell_type] = np.load(a1_scores)
    # load a2 scores
    a2_scores = osp.join(*[
        subcluster_gkmexplain_filepath_base,
        cell_type,
        f"snps.hg38_information.{seq_len}.a2.original.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a2_scores):
        print(f"[WARNING] (Original, subcluster) {cell_type} A2 allele gkmexplain scores does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__original_201bp_gkmexplain_scores_a2_allele[cell_type] = np.load(a2_scores)

''' Shuffled sequence allele specific 201bp gkmexplain scores '''
subcluster__shuffled_201bp_gkmexplain_scores_a1_allele = {}
subcluster__shuffled_201bp_gkmexplain_scores_a2_allele = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    # load a1 scores 
    a1_scores = osp.join(*[
        subcluster_gkmexplain_filepath_base,
        cell_type,
        f"snps.hg38_information.{seq_len}.a1.shuffled.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a1_scores):
        print(f"[WARNING] (Shuffled, subcluster) {cell_type} A1 allele gkmexplain scores does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__shuffled_201bp_gkmexplain_scores_a1_allele[cell_type] = np.load(a1_scores)
    # load a2 scores
    a2_scores = osp.join(*[
        subcluster_gkmexplain_filepath_base,
        cell_type,
        f"snps.hg38_information.{seq_len}.a2.shuffled.imp_score.gkmexplain.normalized.aggregated.npy"
    ])
    if not osp.exists(a2_scores):
        print(f"[WARNING] (Shuffled, subcluster) {cell_type} A2 allele gkmexplain scores does not exist")
        subcluster_cell_types_to_omit += [cell_type]
    else:
        subcluster__shuffled_201bp_gkmexplain_scores_a2_allele[cell_type] = np.load(a2_scores)

''' Original sequence deltaSVM scores '''
subcluster__original_deltaSVM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        subcluster_deltaSVM_filepath_base,
        cell_type, 
        f"snps.hg38_information.{seq_len}.original.{effect_allele}_{non_effect_allele}.deltasvm.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Original, subcluster) {} deltaSVM scores do not exist".format(cell_type))
        subcluster_cell_types_to_omit += [cell_type]
        continue
    subcluster__original_deltaSVM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_deltasvm_score"].tolist())

''' Shuffled sequence deltaSVM scores '''
subcluster__shuffled_deltaSVM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        subcluster_deltaSVM_filepath_base,
        cell_type, 
        f"snps.hg38_information.{seq_len}.shuffled.{effect_allele}_{non_effect_allele}.deltasvm.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Shuffled, subcluster) {} deltaSVM scores do not exist".format(cell_type))
        subcluster_cell_types_to_omit += [cell_type]
        continue
    subcluster__shuffled_deltaSVM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_deltasvm_score"].tolist())

''' Original sequence ISM scores '''
subcluster__original_ISM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        subcluster_ISM_filepath_base,
        cell_type, 
        f"snps.hg38_information.{seq_len}.original.{effect_allele}_{non_effect_allele}.ism.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Original, subcluster) {} ISM scores do not exist".format(cell_type))
        subcluster_cell_types_to_omit += [cell_type]
        continue
    subcluster__original_ISM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_ism_score"].tolist())

''' Shuffled sequence ISM scores '''
subcluster__shuffled_ISM_variant_effect_scores = {}
for broad_cluster, cell_type in zip(broad_clusters, cell_types):
    table =  osp.join(*[
        subcluster_ISM_filepath_base,
        cell_type, 
        f"snps.hg38_information.{seq_len}.shuffled.{effect_allele}_{non_effect_allele}.ism.{cell_type}.aggregated.tsv"
    ])
    if not osp.exists(table):
        print("[WARNING] (Shuffled, subcluster) {} ISM scores do not exist".format(cell_type))
        subcluster_cell_types_to_omit += [cell_type]
        continue
    subcluster__shuffled_ISM_variant_effect_scores[cell_type] = np.array(pd.read_csv(table, sep="\t")["aggregate_ism_score"].tolist())

subcluster_cell_types_to_omit = list(set(subcluster_cell_types_to_omit))

subcluster_data = {
	"original_snp_tables": subcluster__original_snp_tables, # done
	"shuffled_snp_tables": subcluster__shuffled_snp_tables, # done
	"original_201bp_a1_allele_sequence": subcluster__original_201bp_snp_sequences_a1_allele, # done
	"original_201bp_a2_allele_sequence": subcluster__original_201bp_snp_sequences_a2_allele, # done
	"shuffled_201bp_a1_allele_sequence": subcluster__shuffled_201bp_snp_sequences_a1_allele, # done
	"shuffled_201bp_a2_allele_sequence": subcluster__shuffled_201bp_snp_sequences_a2_allele, # done
	"original_gkmexplain_ve_scores": subcluster__original_gkmexplain_variant_effect_scores, # done
	"shuffled_gkmexplain_ve_scores": subcluster__shuffled_gkmexplain_variant_effect_scores, # done
	"original_201bp_gkmexplain_a1_allele_scores": subcluster__original_201bp_gkmexplain_scores_a1_allele, # done
	"original_201bp_gkmexplain_a2_allele_scores": subcluster__original_201bp_gkmexplain_scores_a2_allele, # done
	"shuffled_201bp_gkmexplain_a1_allele_scores": subcluster__shuffled_201bp_gkmexplain_scores_a1_allele, # done
	"shuffled_201bp_gkmexplain_a2_allele_scores": subcluster__shuffled_201bp_gkmexplain_scores_a2_allele, # done
	"original_deltaSVM_ve_scores": subcluster__original_deltaSVM_variant_effect_scores, # done
	"shuffled_deltaSVM_ve_scores": subcluster__shuffled_deltaSVM_variant_effect_scores, # done
	"original_ISM_ve_scores": subcluster__original_ISM_variant_effect_scores, # done
	"shuffled_ISM_ve_scores": subcluster__shuffled_ISM_variant_effect_scores, # done
    "cell_types_to_omit": subcluster_cell_types_to_omit # done
}

''' Sanity check '''
subcluster__cell_types_to_retain = list(set(cell_types) - set(subcluster_data['cell_types_to_omit']))
print("[INFO] The # of cell types to retain for subcluster data is {}".format(len(subcluster__cell_types_to_retain)))
for cell_type in subcluster__cell_types_to_retain:

    # check if these cell types are available within the data structure
    original_data = []
    shuffled_data = []
    for outer_key, data in subcluster_data.items():
        if outer_key == "cell_types_to_omit":
            pass
        elif cell_type not in data:
            print(f"[ERROR, subcluster] {cell_type} is not available in {outer_key} data")
        elif data[cell_type] is None:
            print(f"[ERROR, subcluster] {cell_type} is available in {outer_key} data but empty!")
        elif "original" in outer_key:
            original_data += [data[cell_type]]
        elif "shuffled" in outer_key:
            shuffled_data += [data[cell_type]]
        elif "cell_types_to_omit" in outer_key:
            pass
        else:
            pdb.set_trace()
            print("[ERROR, subcluster] Sanity check (edge case)")

    # check if dimensions of the original data are coherent
    original_data_shapes = [len(x) for x in original_data]
    if len(set(original_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, subcluster] Sanity check (dimensions of the original data are not coherent)")
    else:
        print(f"[INFO, subcluster] {cell_type} original data dimension is {original_data_shapes[0]}")
    
    # check if dimensions of the shuffled data are coherent
    shuffled_data_shapes = [len(x) for x in shuffled_data]
    if len(set(shuffled_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, subcluster] Sanity check (dimensions of the shuffled data are not coherent)")
    else:
        print(f"[INFO, subcluster] {cell_type} shuffled data dimension is {shuffled_data_shapes[0]}")

"""
Save broad and subcluster data
"""
with open(f"{temp_data_dir}/broad_cluster_data.v1.pkl", "wb") as f:
    pickle.dump(broad_cluster_data, f)
with open(f"{temp_data_dir}/subcluster_data.v1.pkl", "wb") as f:
    pickle.dump(subcluster_data, f)

"""
The following part of the script merges broad abd subcluster data
"""

def append_data(data, other_data, other_data_indices):
    """ This function first subsets `other_data` using `other_data_indices` and 
    appends the resulting subset to data for the following data types
        - list
        - np.ndarray
        - pd.DataFrame
    """
    if type(data) == list and type(other_data) == list:
        subset = [other_data[i] for i in other_data_indices]
        return data + subset
    elif type(data) == np.ndarray and type(other_data) == np.ndarray:
        subset = other_data[other_data_indices]
        return np.concatenate((data, subset), axis=0)
    elif type(data) == pd.DataFrame and type(other_data) == pd.DataFrame:
        subset = other_data.iloc[other_data_indices]
        return pd.concat([data, subset], axis=0)
    else:
        raise ValueError("Unsupported data structure type: {}".format(type(data)))

"""
Load broad cluster and subcluster data for each cell type 
"""
with open(f"{temp_data_dir}/broad_cluster_data.v2.pkl", "rb") as f:
    broad_cluster_data = pickle.load(f)
with open(f"{temp_data_dir}/subcluster_data.v2.pkl", "rb") as f:
    subcluster_data = pickle.load(f)

"""
Perform sanity check #1 and #2 to previous step outputs 
"""
# subcluster
# sanity check (1)
for cell_type in subcluster_data.keys():
    # record original and shuffled data separately
    original_data = []
    shuffled_data = []
    for key, data in subcluster_data[cell_type].items():
        if "original" in key:
            original_data.append(data)
        elif "shuffled" in key:
            shuffled_data.append(data)
        else:
            pdb.set_trace()
            print("[ERROR, subcluster] Should not be here")

    # check if the dimensions of the original data are the same
    original_data_shapes = [len(x) for x in original_data]
    if len(set(original_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, subcluster] Sanity check (dimensions of the original data are not coherent)")
    else:
        print(f"[INFO, subcluster] {cell_type} original data dimension is {original_data_shapes[0]}")

    # check if the dimensions of the shuffled data are the same
    shuffled_data_shapes = [len(x) for x in shuffled_data]
    if len(set(shuffled_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, subcluster] Sanity check (dimensions of the shuffled data are not coherent)")
    else:
        print(f"[INFO, subcluster] {cell_type} shuffled data dimension is {shuffled_data_shapes[0]}")

    # check if original_data_dim * 10 == shuffled_data_dim
    if original_data_shapes[0] * null_sequences_count != shuffled_data_shapes[0]:
        pdb.set_trace()
        print("[ERROR, subcluster] Sanity check (original_data_dim * null_sequences_count != shuffled_data_dim)")
    else:
        print(f"[INFO, subcluster] {cell_type} original_data_dim * null_sequences_count == shuffled_data_dim")

# broad cluster 
# Sanity check (2)
for cell_type in broad_cluster_data.keys():
    # record original and shuffled data separately
    original_data = []
    shuffled_data = []
    for key, data in broad_cluster_data[cell_type].items():
        if "original" in key:
            original_data.append(data)
        elif "shuffled" in key:
            shuffled_data.append(data)
        else:
            pdb.set_trace()
            print("[ERROR, broad cluster] Should not be here")

    # check if the dimensions of the original data are the same
    original_data_shapes = [len(x) for x in original_data]
    if len(set(original_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, broad cluster] Sanity check (dimensions of the original data are not coherent)")
    else:
        print(f"[INFO, broad cluster] {cell_type} original data dimension is {original_data_shapes[0]}")

    # check if the dimensions of the shuffled data are the same
    shuffled_data_shapes = [len(x) for x in shuffled_data]
    if len(set(shuffled_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, broad cluster] Sanity check (dimensions of the shuffled data are not coherent)")
    else:
        print(f"[INFO, broad cluster] {cell_type} shuffled data dimension is {shuffled_data_shapes[0]}")

    # check if original_data_dim * 10 == shuffled_data_dim
    if original_data_shapes[0] * null_sequences_count != shuffled_data_shapes[0]:
        pdb.set_trace()
        print("[ERROR, broad cluster] Sanity check (original_data_dim * null_sequences_count != shuffled_data_dim)")
    else:
        print(f"[INFO, broad cluster] {cell_type} original_data_dim * null_sequences_count == shuffled_data_dim")

print("[INFO] Sanity Check (1) and (2) are passed")

# secondary sanity check
for cell_type in broad_cluster_data.keys():
    # broad cluster sanity check 
    broad_original_rsids = broad_cluster_data[cell_type]["original_snp_tables"]["Name"].tolist()
    if len(broad_original_rsids) != len(set(broad_original_rsids)):
        pdb.set_trace()
        print("[ERROR, broad cluster] Sanity check (rsids are not unique)")

    if cell_type in subcluster_data:
        subcluster_original_rsids = subcluster_data[cell_type]["original_snp_tables"]["Name"].tolist()
        if len(subcluster_original_rsids) != len(set(subcluster_original_rsids)):
            pdb.set_trace()
            print("[ERROR, subcluster] Sanity check (rsids are not unique)")

print("Sanity check for broad and subcluster data structures are passed")
        

"""
Merge broad cluster and subcluster data for each cell type
NOTE: broad cluster cell types are a superset of subcluster cell types
Merging strategy 
    - first, set merged original data structures to broad cluster data structures 
    - then, find subcluster snps that are not included in broad cluster structures
    - finally, concatenate subcluster but not broad cluster snps to merged data structures
"""
print("[INFO] merging broad cluster and subcluster data for each cell type")

merged_data = {}
for cell_type in broad_cluster_data.keys():
    # set merged data to broad cluster data
    merged_data[cell_type] = deepcopy(broad_cluster_data[cell_type])

    if cell_type not in subcluster_data:
        continue

    # get original and shuffled snp tables for subcluster and broad clusters
    subcluster_original_snp_table = subcluster_data[cell_type]["original_snp_tables"]
    subcluster_shuffled_snp_table = subcluster_data[cell_type]["shuffled_snp_tables"]
    broad_cluster_original_snp_table = broad_cluster_data[cell_type]["original_snp_tables"]
    broad_cluster_shuffled_snp_table = broad_cluster_data[cell_type]["shuffled_snp_tables"]

    # get snps from the broad and subcluster original tables
    subcluster_snps = subcluster_original_snp_table["Name"].tolist()
    broad_cluster_snps = broad_cluster_original_snp_table["Name"].tolist()

    # perform sanity check for uniqueness of the snps in broad and subclusters
    if len(broad_cluster_snps) != len(set(broad_cluster_snps)):
        pdb.set_trace()
        print("[ERROR, broad cluster] Sanity check (rsids are not unique)")
    if len(subcluster_snps) != len(set(subcluster_snps)):
        pdb.set_trace()
        print("[ERROR, subcluster] Sanity check (rsids are not unique)")
    

    # find snps that are in subcluster but not broad clusters 
    # and check if importance scores of these SNPs are the same
    # in broad and subcluster data structures 
    sub_diff_broad_original_indices = []
    for subcluster_idx, rsid in enumerate(subcluster_snps):
        try:
            # find snps in subcluster that are also in broad cluster
            broad_idx = broad_cluster_snps.index(rsid)
            # compare scores between broad and subcluster for the overlapping snps 
            assert subcluster_data[cell_type]["original_gkmexplain_ve_scores"][subcluster_idx] == \
                broad_cluster_data[cell_type]["original_gkmexplain_ve_scores"][broad_idx]
            assert subcluster_data[cell_type]["original_deltaSVM_ve_scores"][subcluster_idx] == \
                broad_cluster_data[cell_type]["original_deltaSVM_ve_scores"][broad_idx]
            assert subcluster_data[cell_type]["original_ISM_ve_scores"][subcluster_idx] == \
                broad_cluster_data[cell_type]["original_ISM_ve_scores"][broad_idx]
        except ValueError:
            # the rsid is not available in broad clusters, record it for
            # appending the data to merged data structures 
            sub_diff_broad_original_indices.append(subcluster_idx)
        except AssertionError:
            pdb.set_trace()
            print("Scores of the snp does not match between broad and subcluster")

    # extend the indices for original data to broad data
    sub_diff_broad_shuffled_indices = []
    for idx in sub_diff_broad_original_indices:
        for i in range(null_sequences_count):
            sub_diff_broad_shuffled_indices += [idx * null_sequences_count + i]
    
    # concatenate sub diff broad snps to merged data structures 
    for key in merged_data[cell_type]:
        if "shuffled" in key:
            merged_data[cell_type][key] = \
                append_data(merged_data[cell_type][key], subcluster_data[cell_type][key], sub_diff_broad_shuffled_indices)
        elif "original" in key:
            merged_data[cell_type][key] = \
                append_data(merged_data[cell_type][key], subcluster_data[cell_type][key], sub_diff_broad_original_indices)
        else:
            pdb.set_trace()
            print("Should not be here")


"""
Sanity check #3 to the merged data structures
"""
for cell_type in merged_data.keys():
    # record original and shuffled data separately
    original_data = []
    shuffled_data = []
    for key, data in merged_data[cell_type].items():
        if "original" in key:
            original_data.append(data)
        elif "shuffled" in key:
            shuffled_data.append(data)
        else:
            pdb.set_trace()
            print("[ERROR, merged] Should not be here")

    # check if the dimensions of the original data are the same
    original_data_shapes = [len(x) for x in original_data]
    if len(set(original_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, merged] Sanity check (dimensions of the original data are not coherent)")
    else:
        print(f"[INFO, merged] {cell_type} original data dimension is {original_data_shapes[0]}")

    # check if the dimensions of the shuffled data are the same
    shuffled_data_shapes = [len(x) for x in shuffled_data]
    if len(set(shuffled_data_shapes)) > 1:
        pdb.set_trace()
        print("[ERROR, merged] Sanity check (dimensions of the shuffled data are not coherent)")
    else:
        print(f"[INFO, merged] {cell_type} shuffled data dimension is {shuffled_data_shapes[0]}")

    # (3) check if original_data_dim * 10 == shuffled_data_dim
    if original_data_shapes[0] * null_sequences_count != shuffled_data_shapes[0]:
        pdb.set_trace()
        print("[ERROR, merged] Sanity check (original_data_dim * null_sequences_count != shuffled_data_dim)")
    else:
        print(f"[INFO, merged] {cell_type} original_data_dim * null_sequences_count == shuffled_data_dim")

print("[INFO] Sanity check #3 is passed")


'''
Sanity check #2 and #3 to the merged data structures
'''
for cell_type in merged_data.keys():
    original_snp_table = merged_data[cell_type]["original_snp_tables"]
    shuffled_snp_table = merged_data[cell_type]["shuffled_snp_tables"]

    original_rsids = original_snp_table["Name"].tolist()
    shuffled_rsids = shuffled_snp_table["Name"].tolist()

    # (4) check if all snps in original rsids are in shuffled snp list
    if len(original_rsids) != len(set(original_rsids)):
        pdb.set_trace()
        print("[ERROR, merged] Sanity check (original rsids are not unique)")

    
    
    for idx, rsid in enumerate(original_rsids):
        bool_vector = [True if rsid == x else False for x in shuffled_rsids]
        # (5.1) check if number of occurences of original rsids in shuffled snp list equals to null_sequences_count
        assert sum(bool_vector) == null_sequences_count
        # (5.2) check if all occurences are at their expected place
        expected_match_locs = range(idx * null_sequences_count, (idx + 1) * null_sequences_count)
        bool_vector_subset = [bool_vector[idx] for idx in expected_match_locs]
        assert sum(bool_vector_subset) == null_sequences_count

print("[INFO] Sanity check #4 and #5 are passed")

"""
Save merged data structure
"""
with open(f"{temp_data_dir}/merged_data.v1.pkl", "wb") as f:
    pickle.dump(merged_data, f)


print("Script finished")





print("Script finished")
