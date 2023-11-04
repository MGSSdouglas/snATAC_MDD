import os 
import pdb
import sys
import pickle
import random
import argparse 

import numpy as np 
import os.path as osp
import pandas as pd 
import numpy as np 
import pyranges as pr

from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO


seed = 66
random.seed(seed)
np.random.seed(seed)


# ------- Command line argument parsing -------

parser = argparse.ArgumentParser(description="Generate inaccessible sequences")

parser.add_argument("-c", "--cluster", dest="cluster", type=str,
                    help="Cluster")
parser.add_argument("-o", "--ocr-dir", dest="ocr_dir", type=str,
                    help="Directory of cluster-specific OCRs (.tsv)")
parser.add_argument("-n", "--candidate-negative-file", dest="cand_filepath", type=str,
                    help="Candidate negative genome regions file (.tsv)")
parser.add_argument("-t", "--target-dir", dest="target_dir", type=str,
                    help="Target directory where dataset will be saved")
parser.add_argument("-g", "--genome-path", dest="genome_filepath", type=str,
                    help="Path to genome dictionary file (.pkl)")
parser.add_argument("-fs", "--fold-split-file", dest="fold_splits", type=str,
                    help="Chromosome to CV fold mapping file (.txt)")
parser.add_argument("-ngcb", "--num-gc-bins", dest="num_gc_bins", type=int,
                    help="Number of equally populated GC bins")

args = parser.parse_args()

# ------- Utils -------

def find_gc_bin_membership(gc, gc_bin_edges):
    count = 0
    for i in range(len(gc_bin_edges)-1):
        if i == 0:
            if gc < gc_bin_edges[i]:
                count = -1
                break
            elif gc < gc_bin_edges[i+1]:
                break
            else:
                count += 1
        elif i != len(gc_bin_edges) - 2:
            if gc >= gc_bin_edges[i] and gc < gc_bin_edges[i+1]:
                break
            else:
                count += 1
        else:
            if gc >= gc_bin_edges[i] and gc <= gc_bin_edges[i+1]:
                break
            else:
                count = -1
                break
    return count

# load prepared hg38 reference genome
with open(args.genome_path, "rb") as f:
    genome = pickle.load(f)

extract_seq = lambda chr, start, end: genome[chrom][start:end+1]
compute_GC = lambda sequence: GC(sequence) / 100
gc_bin_membership = lambda x: find_gc_bin_membership(x, gc_bin_edges)
get_record = lambda x: SeqRecord(Seq(x["sequence"]), id=x["header"], description="")

# ------- Script -------

# load OCRs for the cluster
filepath = osp.join(args.ocr_dir, f"{args.cluster}_OCR.tsv")
ocrs = pd.read_csv(filepath, sep="\t")


# extract underlying ref genome sequence for OCRs
ocrs["Sequence"] = ocrs[["Chromosome", "Start", "End"]].map(extract_seq)

# calculate GC content for OCRs
ocrs["GC_1001bp"] = np.around(ocrs["Sequence"].map(compute_GC), decimals=4)

# OCR GC content bins
num_gc_bins = args.num_gc_bins
step_size = int(100/num_gc_bins)
percentiles = np.arange(0, 100+step_size, step_size)
gc_bin_edges = np.percentile(ocrs["GC_1001bp"], percentiles)

# OCR GC bin membership
ocrs["GC_bin"] = ocrs["GC_1001bp"].map(gc_bin_membership)

# retain relevant columns
ocrs = ocrs[["neglog_MACS2_pval", "GC_bin", "Chromosome", "Start", "End", "peak_id", "absolute_peak_id", "Sequence", "Strand"]]

# convert to pyranges
pr_ocrs = pr.PyRanges(ocrs, int64=True)

# load candidate negative genome regions
candidate_negs = pd.read_csv(args.cand_filepath, sep="\t")
candidate_negs = candidate_negs[["Chromosome", "Start", "End", "GC_1001bp"]]
candidate_negs.rename(columns={"GC_1001bp": "cand_GC_1001bp"}, inplace=True)

# convert to pyranges
pr_candidate_negs = pr.PyRanges(candidate_negs, int64=True)

# determine coverage of candidate negative regions on open chromatin regions
overlaps = pr_candidate_negs.coverage(pr_ocrs, strandedness=False,
                overlap_col="C", fraction_col="F", nb_cpu=4).as_df()

# retain candidate negative peaks with no overlap only
overlaps = overlaps[overlaps["F"] == 0]
overlaps["Strand"] = ["+"] * overlaps.shape[0]
overlaps["distance_u"] = [-1] * overlaps.shape[0]
overlaps["distance_d"] = [-1] * overlaps.shape[0]
overlaps["Index"] = list(range(overlaps.shape[0]))

overlaps["cand_GC_bin"] = overlaps["cand_GC_1001bp"].map(gc_bin_membership)

overlaps = overlaps[["Chromosome", "Start", "End", "cand_GC_bin"]]

# remove peaks with GC bin membership -1 (not in the range of positive peaks)
overlaps = overlaps[overlaps["cand_GC_bin"] != -1]


# ------- Prepraration before inaccessible sequence sampling -------

# split negative samples with respect to chromosomes
# single level dictionary: {chromosome: pd.DataFrame of peaks}
gps = overlaps.groupby("Chromosome")
cand_negative_dct = {chrom:gps.get_group(chrom).copy() for chrom in gps.groups}
chroms = list(cand_negative_dct.keys())
del gps, overlaps

# three level dictionary: {chromosome: {gc_bin_membership: {"population":pd.DataFrame of peaks, "next_sample_idx":int}}}
temp = {chrom:{} for chrom in chroms}
# iterate over chromosomes
for chrom in chroms:
    table = cand_negative_dct[chrom]
    gps = table.groupby("GC_bin")
    # iterate over GC bins
    for gc_bin in gps.groups:
        table_ = gps.get_group(gc_bin).copy()
        # shuffle candidate negative samples
        table_ = table_.iloc[np.random.permutation(len(table_))]
        # reflect changes
        temp[chrom][gc_bin] = {}
        temp[chrom][gc_bin]["population"] = table_
        temp[chrom][gc_bin]["next_sample_idx"] = 0
    del cand_negative_dct[chrom]
del cand_negative_dct


# ------- Negative (inaccessible) sequence sampling per OCR -------

cand_negative_dct = temp
negative_dct = dict(Chromosome=[], Start=[], End=[], GC_bin=[], peak_id=[])

# iterate over positive examples (OCRs)
for index, row in ocrs.iterrows():

    # select the negative peak set which matches chromosome and GC bin with current positive peak
    temp = cand_negative_dct[row["Chromosome"]][row["GC_bin"]]
    table = temp["population"]
    next_sample_idx = temp["next_sample_idx"]

    if next_sample_idx == len(table):
        raise ValueError(f"No more negative peaks to sample for GC bin {row['GC_bin']} !")
    sampled_peak = table.iloc[next_sample_idx]

    negative_dct["Chromosome"].append(sampled_peak["Chromosome"])
    negative_dct["Start"].append(sampled_peak["Start"])
    negative_dct["End"].append(sampled_peak["End"])
    negative_dct["peak_id"].append(row["peak_id"])
    negative_dct["GC_bin"].append(row["GC_bin"])

    # update last sampled index
    cand_negative_dct[row["Chromosome"]][row["GC_bin"]]["next_sample_idx"] += 1

# convert to pd.DataFrame
negatives = pd.DataFrame(negative_dct)


# ------- Combine positives and negatives -------

positives = ocrs

# add labels
negatives["label"] = 0
positives["label"] = 1

# set indices
negatives.set_index("peak_id", inplace=True)
positives.set_index("peak_id", inplace=True)

# retain relevant cols
positives = positives[["label", "GC_bin", "Chromosome", "Start", "End", "neglog_MACS2_pval"]]
negatives = negatives[["label", "GC_bin", "Chromosome", "Start", "End"]]

# concate positives and negatives
negatives["neglog_MACS2_pval"] = -1 
dataset = pd.concat([positives, negatives], axis=0)


# ------- Partition positive and negative samples to cross validation folds -------

# read cv fold chromosome contents
with open(args.fold_split_file, "r") as f:
    fold2chroms = {f"fold_{idx+1}": set(val.strip().split(",")) \
        for idx, val in enumerate(f.readlines()) if val.strip() != ''}
    chromosomes = set.union(*fold2chroms.values())   

# generate two-level dataset dct: fold_id -> train/test 
cv_dataset = {}
for fold_id, fold_chroms in fold2chroms.items():
    other_chroms = chromosomes - fold_chroms
    cv_dataset[fold_id] = {
        "test": dataset[dataset.Chromosome.str.isin(fold_chroms)],
        "train": dataset[dataset.Chromosome.str.isin(other_chroms)]
    }

# ------- Retain top 60,000 pos-neg training example pairs and recycle unused candidate negs ------- 

rank_thr = 60000
for fold_id in cv_dataset:

    train_data = cv_dataset[fold_id]["train"]

    train_data_positives = train_data[train_data.neglog_MACS2_pval != -1]

    train_data_positives.sort_values(by="neglog_MACS2_pval", ascending=False)

    # retain top 60K (positive, negative) training example pairs wrt MACS2 p-value
    top_60K = train_data_positives.iloc[:rank_thr,:].index
    cv_dataset[fold_id]["train"] = \
        train_data[train_data.index.str.isin(top_60K)]

    # extract unused negative examples
    rest = train_data_positives.iloc[rank_thr:,:].index
    rest_negatives = train_data[
        (train_data.index.str.isin(rest)) &
        (train_data.label == 0)
    ]

    # recycle extracted examples to candidate negative examples
    for chrom in cand_negative_dct:
        for gc_bin in cand_negative_dct[chrom]:

            recycled_negatives = rest_negatives[
                (rest_negatives.Chromosome == chrom) &
                (rest_negatives.GC_bin == gc_bin)
            ]
            recycled_negatives = set(recycled_negatives.peak_id.tolist())

            temp = cand_negative_dct[row["Chromosome"]][row["GC_bin"]]
        
            population = set(temp["population"].peak_id.tolist())
            next_sample_idx = temp["next_sample_idx"]

            right = population.iloc[next_sample_idx:,:]
            right = set(right.peak_id_)

            left = population.iloc[:next_sample_idx,:]
            left = set(left.peak_id_)

            try:
                assert len(recycled_negatives & left) == len(recycled_negatives)
            except:
                pdb.set_trace()
            
            temp = temp[temp.peak_id.str.isin(recycled_negatives+right)]
            
            cand_negative_dct[row["Chromosome"]][row["GC_bin"]]["population"] = \
                temp.iloc[np.random.permutation(len(temp))]
            cand_negative_dct[row["Chromosome"]][row["GC_bin"]]["next_sample_idx"] = 0
            
    
# ------- Split pos, neg and generate fasta ------- 


for fold_id in cv_dataset:

    target_dir = osp.join(args.target_dir, f"{args.cluster}/{fold_id}")
    os.makedirs(target_dir, exist_ok=True)

    # train split

    train_data = cv_dataset[fold_id]["train"]
    filepath = osp.join(target_dir, "CV.train.tsv")
    train_data.to_csv(filepath, sep="\t")

    train_positives = train_data[train_data.label == 1]
    filepath = osp.join(target_dir, "CV.train.pos.tsv")
    train_positives.to_csv(filepath, sep="\t")

    train_negatives = train_data[train_data.label == 0]
    filepath = osp.join(target_dir, "CV.train.neg.tsv")
    train_negatives.to_csv(filepath, sep="\t")

    train_positives["Sequence"] = \
        train_positives[["Chromosome", "Start", "End"]].map(extract_seq).str.upper()
    
    train_negatives["Sequence"] = \
        train_negatives[["Chromosome", "Start", "End"]].map(extract_seq).str.upper()
    
    filepath = osp.join(target_dir, "CV.train.positives.fa")
    with open(filepath, "w") as f:

        # https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        f_out = FastaIO.FastaWriter(f, wrap=None)

        records = [get_record({"header": peak_id, "sequence": seq})
            for peak_id, seq in zip(train_positives["peak_id"], train_positives["Sequence"])]
        f_out.write_file(records)

    filepath = osp.join(target_dir, "CV.train.negatives.fa")
    with open(filepath, "w") as f:

        # https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        f_out = FastaIO.FastaWriter(f, wrap=None)

        records = [get_record({"header": peak_id, "sequence": seq})
            for peak_id, seq in zip(train_negatives["peak_id"], train_negatives["Sequence"])]
        f_out.write_file(records)


    # test split

    test_data = cv_dataset[fold_id]["train"]
    filepath = osp.join(target_dir, "CV.test.tsv")
    test_data.to_csv(filepath, sep="\t")

    test_positives = test_data[test_data.label == 1]
    filepath = osp.join(target_dir, "CV.test.pos.tsv")
    test_positives.to_csv(filepath, sep="\t")

    test_negatives = test_data[test_data.label == 0]
    filepath = osp.join(target_dir, "CV.test.neg.tsv")
    test_negatives.to_csv(filepath, sep="\t")

    test_positives["Sequence"] = \
        test_positives[["Chromosome", "Start", "End"]].map(extract_seq).str.upper()
    
    test_negatives["Sequence"] = \
        test_negatives[["Chromosome", "Start", "End"]].map(extract_seq).str.upper()
    
    filepath = osp.join(target_dir, "CV.test.pos.fa")
    with open(filepath, "w") as f:

        # https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        f_out = FastaIO.FastaWriter(f, wrap=None)

        records = [get_record({"header": peak_id, "sequence": seq})
            for peak_id, seq in zip(test_positives["peak_id"], test_positives["Sequence"])]
        f_out.write_file(records)

    filepath = osp.join(target_dir, "CV.test.neg.fa")
    with open(filepath, "w") as f:

        # https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        f_out = FastaIO.FastaWriter(f, wrap=None)

        records = [get_record({"header": peak_id, "sequence": seq})
            for peak_id, seq in zip(test_negatives["peak_id"], test_negatives["Sequence"])]
        f_out.write_file(records)


# ------- Negative (inaccessible) sequence sampling per marker cCRE -------

# load marker cCREs 
filepath = f"{args.ocr_dir}/../cluster__marker_cCREs.tsv"
marker_ccres = pd.read_csv(filepath, sep="\t")

# retain relevant columns
marker_ccres = marker_ccres[["Chromosome", "Start", "End", "cluster_id", "peak_id"]]

# extract sequence -> calculate GC -> GC bin membership
marker_ccres["GC_bin"] = \
    marker_ccres[["Chromosome", "Start", "End"]].map(extract_seq).map(compute_GC).map(gc_bin_membership)

# iterate over positive examples (marker cCREs)
negative_dct = dict(Chromosome=[], Start=[], End=[], GC_bin=[], peak_id=[])
for index, row in marker_ccres.iterrows():

    # select the negative peak set which matches chromosome and GC bin with current positive peak
    temp = cand_negative_dct[row["Chromosome"]][row["GC_bin"]]
    table = temp["population"]
    next_sample_idx = temp["next_sample_idx"]

    if next_sample_idx == len(table):
        raise ValueError(f"No more negative sequences to sample for GC bin {row['GC_bin']} !")
    sampled_peak = table.iloc[next_sample_idx]

    negative_dct["Chromosome"].append(sampled_peak["Chromosome"])
    negative_dct["Start"].append(sampled_peak["Start"])
    negative_dct["End"].append(sampled_peak["End"])
    negative_dct["peak_id"].append(row["peak_id"])
    negative_dct["GC_bin"].append(row["GC_bin"])

    # update last sampled index
    cand_negative_dct[row["Chromosome"]][row["GC_bin"]]["next_sample_idx"] += 1

# convert to pd.DataFrame
negatives = pd.DataFrame(negative_dct)


# ------- (marker cCRE) Combine positive and negative examples -------

positives = marker_ccres

# add labels
negatives["label"] = 0
positives["label"] = 1

# set indices
negatives.set_index("peak_id", inplace=True)
positives.set_index("peak_id", inplace=True)

# retain relevant cols
positives = positives[["label", "GC_bin", "Chromosome", "Start", "End"]]
negatives = negatives[["label", "GC_bin", "Chromosome", "Start", "End"]]
dataset = pd.concat([positives, negatives], axis=0)


# ------- (marker cCRE) Partition positive and negative examples to cross validation folds wrt chromosomes-------

# read cv fold chromosome contents
with open(args.fold_split_file, "r") as f:
    fold2chroms = {f"fold_{idx+1}": set(val.strip().split(",")) \
        for idx, val in enumerate(f.readlines()) if val.strip() != ''}
    chromosomes = set.union(*fold2chroms.values())   

# split marker cCRE dataset according to fold specifications
marker_eval_dataset = {
    fold_id: marker_eval_dataset[marker_eval_dataset.Chromosome.str.isin(fold_chroms)]
        for fold_id, fold_chroms in fold2chroms.items()
}


# ------- (marker cCRE) Split pos, neg and generate fasta -------

for fold_id in cv_dataset:

    target_dir = osp.join(args.target_dir, f"{args.cluster}/{fold_id}")
    os.makedirs(target_dir, exist_ok=True)

    # train split

    data = cv_dataset[fold_id]["train"]
    filepath = osp.join(target_dir, "marker.test.tsv")
    data.to_csv(filepath, sep="\t")

    positives = data[data.label == 1]
    filepath = osp.join(target_dir, "marker.test.pos.tsv")
    positives.to_csv(filepath, sep="\t")

    negatives = data[data.label == 0]
    filepath = osp.join(target_dir, "marker.test.neg.tsv")
    negatives.to_csv(filepath, sep="\t")

    positives["Sequence"] = \
        positives[["Chromosome", "Start", "End"]].map(extract_seq).str.upper()
    
    negatives["Sequence"] = \
        negatives[["Chromosome", "Start", "End"]].map(extract_seq).str.upper()
    
    filepath = osp.join(target_dir, "marker.test.pos.fa")
    with open(filepath, "w") as f:

        # https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        f_out = FastaIO.FastaWriter(f, wrap=None)

        records = [get_record({"header": peak_id, "sequence": seq})
            for peak_id, seq in zip(positives["peak_id"], positives["Sequence"])]
        f_out.write_file(records)

    filepath = osp.join(target_dir, "marker.test.neg.fa")
    with open(filepath, "w") as f:

        # https://stackoverflow.com/questions/24156578/using-bio-seqio-to-write-single-line-fasta
        f_out = FastaIO.FastaWriter(f, wrap=None)

        records = [get_record({"header": peak_id, "sequence": seq})
            for peak_id, seq in zip(negatives["peak_id"], negatives["Sequence"])]
        f_out.write_file(records)

print("Script finished")
