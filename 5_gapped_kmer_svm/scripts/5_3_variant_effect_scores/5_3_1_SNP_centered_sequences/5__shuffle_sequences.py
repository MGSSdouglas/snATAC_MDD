import pdb 
import os 
import pickle
import subprocess

import os.path as osp
import numpy as np 
import pandas as pd 

import pyranges as pr

from Bio import SeqIO
from tqdm import tqdm

'''
This script generates given number of dinucleotide shuffled sequences
from 51bp and 201bp reference genome sequences, centered at cdSN.

NOTE: load the following modules before loading python environment and therefore running this script:
    -   module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
    -   module load meme/5.4.1
'''

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# config for shuffle
num_shuffle = 10
shuffle_suffix = []
random_seed = 35
seqlens = [51, 201]

# get cdSNPs basepath
cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")

# iterate over clusters
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

    # iterate over sequence lengths
    for seqlen in seqlens:

        fpath = osp.join(*[cdsnp_basepath, cluster, f"{seqlen}bp_seqs"])

        no_allele_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.no_allele.original.fa"])
        shuffled_no_allele_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.no_allele.shuffled.fa"])

        '''
        Step 1.a: Dinuc preserved shuffle 10 times and save resulting sequences
        '''
        # command to call fasta-dinucleotide-shuffle from MEME Suite
        cmd = ["fasta-dinucleotide-shuffle"]
        cmd += ["-f", no_allele_fasta_file] # input fasta file
        cmd += ["-s", f"{random_seed}"] # random seed for shuffling
        cmd += ["-c", str(num_shuffle)] # number of shuffles
        # output the resulting fasta to a file 
        cmd += [">", shuffled_no_allele_fasta_file]

        # submit job
        print(' '.join(cmd))
        result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
        if result.returncode != 0:
            pdb.set_trace()

        '''
        Step 1.b: introduce middle position (A1 and A2, separately) back into sequences
        '''

        middle_pos = seqlen // 2

        a1_alleles = snps_in_peaks["A1"].values.tolist()
        a2_alleles = snps_in_peaks["A2"].values.tolist()
        rsids = snps_in_peaks["Name"].values.tolist()

        # [item1, item2] -> [num_shuffle * [item1], num_shuffle * [item2]]
        a1_alleles = [a1_alleles[i//num_shuffle] for i in range(num_shuffle*len(a1_alleles))]
        a2_alleles = [a2_alleles[i//num_shuffle] for i in range(num_shuffle*len(a2_alleles))]
        rsids =  [rsids[i//num_shuffle] for i in range(num_shuffle*len(rsids))]

        shuffled_a1_fasta = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.a1.shuffled.fa"])
        with open(shuffled_a1_fasta, "w") as a1, \
            open(shuffled_no_allele_fasta_file, "r") as no_allele:

            for idx, record in tqdm(enumerate(SeqIO.parse(no_allele, "fasta")), desc="A1", leave="False"):

                # insert a1 allele into the middle position
                record.seq = record.seq[:middle_pos] + a1_alleles[idx] + record.seq[middle_pos:]
                record.description = rsids[idx]
                SeqIO.write(record, a1, "fasta")

        shuffled_a2_fasta = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.a1.shuffled.fa"])
        with open(shuffled_a2_fasta, "w") as a2, \
            open(shuffled_no_allele_fasta_file, "r") as no_allele:

            for idx, record in tqdm(enumerate(SeqIO.parse(no_allele, "fasta")), desc="A2", leave="False"):

                # insert a2 allele into the middle position
                record.seq = record.seq[:middle_pos] + a2_alleles[idx] + record.seq[middle_pos:]
                record.description = rsids[idx]
                SeqIO.write(record, a2, "fasta")

