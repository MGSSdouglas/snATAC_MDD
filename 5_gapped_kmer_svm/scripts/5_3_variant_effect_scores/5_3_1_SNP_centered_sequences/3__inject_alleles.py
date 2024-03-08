import pdb 
import os 
import pickle

import os.path as osp
import numpy as np 
import pandas as pd 

import pyranges as pr

from Bio import SeqIO

"""
This script injects cdSNP alleles into their repspective positions in
the 51bp and 201bp reference genome sequences centered at cdSNPs
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
    for seqlen in seqlens:

        fpath = osp.join(*[cdsnp_basepath, cluster, f"{seqlen}bp_seqs"])

        # find the index of the center position
        middle_pos = seqlen // 2

        # copy reference fasta file to A1, A2 and no_allele fasta files
        ref_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.ref.fa"])
        a1_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.a1.original.fa"])
        a2_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.a2.original.fa"])
        no_allele_fasta_file = osp.join(*[fpath, f"cdsnps.unique.{seqlen}bp.no_allele.original.fa"])

        # open ref and a1 fasta files and mutate the middle location of ref sequences to add a1 alleles
        with open(ref_fasta_file, "r") as ref, open(a1_fasta_file, "w") as a1:
            for idx, record in enumerate(SeqIO.parse(ref, "fasta")):
                record.seq = record.seq[:middle_pos] + a1_alleles[idx] + record.seq[middle_pos + 1:]
                record.description = snp_rsids[idx]
                SeqIO.write(record, a1, "fasta")

        # open ref and a2 fasta files and mutate the middle location of ref sequences to add a2 alleles
        with open(ref_fasta_file, "r") as ref, open(a2_fasta_file, "w") as a2:
            for idx, record in enumerate(SeqIO.parse(ref, "fasta")):
                record.seq = record.seq[:middle_pos] + a2_alleles[idx] + record.seq[middle_pos + 1:]
                record.description = snp_rsids[idx]
                SeqIO.write(record, a2, "fasta")

        # open ref and noallele fasta files and remove the middle location of ref sequences
        with open(ref_fasta_file, "r") as ref, open(no_allele_fasta_file, "w") as no_allele:
            for idx, record in enumerate(SeqIO.parse(ref, "fasta")):
                record.seq = record.seq[:middle_pos] + record.seq[middle_pos + 1:]
                assert len(record.seq) == seqlen - 1
                record.description = snp_rsids[idx]
                SeqIO.write(record, no_allele, "fasta")


