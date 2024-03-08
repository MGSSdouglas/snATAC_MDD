import pdb 
import os 
import pickle

import os.path as osp
import numpy as np 
import pandas as pd 

import pyranges as pr

from Bio import SeqIO

"""
This script perform sanity check by comparing 51bp and 201bp sequences
to see if alleles are injected to the correct positions
"""

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# get cdSNPs basepath
cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")

# iterate over clusters and generate 51bp and 201bp allelic sequences for cdSNPs
seqlens = [51, 201]
clusters = os.listdir(cdsnp_basepath)
for cluster in clusters:

    # load cdsnps for the cluster
    try:
        snps_in_peaks = pd.read_csv( \
            osp.join(*[cdsnp_basepath, cluster, "cdsnps.unique.tsv"]), sep="\t"
        )
    except pd.errors.EmptyDataError:
        print(f"EmptyDataError: {cluster}")
        continue

    # get a1 and a2 alleles
    a1_alleles = snps_in_peaks["A1"]
    a2_alleles = snps_in_peaks["A2"]
    snp_rsids = snps_in_peaks["Name"]

    # iterate over seqlens
    a1_sequences = {}
    a2_sequences = {}
    for seqlen in seqlens:

        fpath = osp.join(*[cdsnp_basepath, cluster, f"{seqlen}bp_seqs"])

        # Extract sequences for A1 allele of cdSNPs
        a1_sequences[f"{seqlen}bp"] = {}
        a1_sequences[f"{seqlen}bp"]["rsids"] = []
        a1_sequences[f"{seqlen}bp"]["sequences"] = []
        a1_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.a1.original.fa"])
        with open(a1_fasta_file, "r") as a1:
            for record in SeqIO.parse(a1, "fasta"):
                a1_sequences[f"{seqlen}bp"]["rsids"].append(record.description)
                a1_sequences[f"{seqlen}bp"]["sequences"].append(record.seq)

        # Extract sequences for A2 allele of cdSNPs
        a2_sequences[f"{seqlen}bp"] = {}
        a2_sequences[f"{seqlen}bp"]["rsids"] = []
        a2_sequences[f"{seqlen}bp"]["sequences"] = []
        a2_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.a2.original.fa"])
        with open(a2_fasta_file, "r") as a2:
            for record in SeqIO.parse(a2, "fasta"):
                a2_sequences[f"{seqlen}bp"]["rsids"].append(record.description)
                a2_sequences[f"{seqlen}bp"]["sequences"].append(record.seq)
        
    #### Test 1: compare middle position of 51bp and 201bp sequences 

    ## A1 allele
    for idx, (snp, a1) in enumerate(zip(snp_rsids, a1_alleles)):
        alleles = []
        for seqlen in seqlens:
            middle_pos = seqlen // 2
            alleles.append(a1_sequences[f"{seqlen}bp"]["sequences"][idx][middle_pos])
            cur_snp = a1_sequences[f"{seqlen}bp"]["rsids"][idx].split(" ")[-1]
            assert cur_snp == snp, f" (A1) {seqlen}bp snp {cur_snp} != {snp}"

        print(f"(A1) {snp}, {a1}, {alleles}")
        
        for seqlen, allele in zip(seqlens,alleles):
            if allele != a1:
                print(f"(A1) {seqlen}bp sequence for {snp} has allele {allele}, but should have been {a1}")

    ## A2 allele

    for idx, (snp, a2) in enumerate(zip(snp_rsids, a2_alleles)):
        alleles = []
        for seqlen in seqlens:
            middle_pos = seqlen // 2
            alleles.append(a2_sequences[f"{seqlen}bp"]["sequences"][idx][middle_pos])
            cur_snp = a2_sequences[f"{seqlen}bp"]["rsids"][idx].split(" ")[-1]
            assert cur_snp == snp, f" (A2) {seqlen}bp snp {cur_snp} != {snp}"

        print(f"(A2) {snp}, {a2}, {alleles}")

        for seqlen, allele in zip(seqlens,alleles):
            if allele != a2:
                print(f"(A2) {seqlen}bp sequence for {snp} has allele {allele}, but should have been {a2}")

    #### Test 2: check if sequences are consistent across different sequence lengths 

    seqlen_pairs = [(51, 201)]

    for pair in seqlen_pairs:
        source_seqlen = pair[0]
        target_seqlen = pair[1]

        source_middle_pos = source_seqlen // 2
        target_middle_pos = target_seqlen // 2

        ## A1
        for idx, (snp, a1) in enumerate(zip(snp_rsids, a1_alleles)):
            source = a1_sequences[f"{source_seqlen}bp"]["sequences"][idx]
            target = a1_sequences[f"{target_seqlen}bp"]["sequences"][idx]
            if source != target[target_middle_pos - source_middle_pos : target_middle_pos + source_middle_pos + 1]:
                pdb.set_trace()
        ## A2
        for idx, (snp, a2) in enumerate(zip(snp_rsids, a2_alleles)):
            source = a2_sequences[f"{source_seqlen}bp"]["sequences"][idx]
            target = a2_sequences[f"{target_seqlen}bp"]["sequences"][idx]
            if source != target[target_middle_pos - source_middle_pos : target_middle_pos + source_middle_pos + 1]:
                pdb.set_trace()

