import pdb 
import os 
import os.path as osp

import numpy as np 
import pandas as pd 

from tqdm import tqdm

'''
This script splits the sequences into chunks and saves them to files.
'''

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# configuration
num_chunks = 16
target_seqlens = [51, 201]
target_snp_types = ["a1", "a2"]
target_seq_types = ["original", "shuffled"]

# get cdSNPs basepath
cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")

# iterate over clusters
clusters = os.listdir(cdsnp_basepath)
for cluster in clusters:
    filepath_base = osp.join(*[cdsnp_basepath, cluster])

    for target_seqlen in tqdm(target_seqlens, desc="seqlens"):
        for target_snp_type in tqdm(target_snp_types, desc="snp types", leave=False):
            for target_seq_type in tqdm(target_seq_types, desc="seq types", leave=False):

                print(f"> Splitting {target_seqlen}bp, {target_snp_type}, {target_seq_type} seqs to {num_chunks} chunks")

                if target_seq_type == "original":
                    input_filename = f"cdsnps.unique.{target_seqlen}bp.{target_snp_type}.original"
                elif target_seq_type == "shuffled":
                    input_filename = f"cdsnps.unique.{target_seqlen}bp.{target_snp_type}.shuffled"

                input_dir = f"{filepath_base}/{target_seqlen}bp_seqs"
                input_filepath = f"{input_dir}/{input_filename}.fa"

                if not osp.exists(input_filepath):
                    print(f"File not found: {input_filepath}")
                    continue

                chunk_dir = f"{filepath_base}/{target_seqlen}bp_seqs/chunks/{target_snp_type}/{target_seq_type}/splits"
                os.makedirs(chunk_dir, exist_ok=True)

                merged_dir = f"{filepath_base}/{target_seqlen}bp_seqs/chunks/{target_snp_type}/{target_seq_type}/merged"
                os.makedirs(merged_dir, exist_ok=True)

                # read input fasta file to a list where each element is a fasta record
                with open(input_filepath, "r") as f:
                    lines = f.readlines()

                records = []
                for line in lines:
                    if line.startswith(">"):
                        records.append(line)
                    else:
                        records[-1] += line

                # split records into {num_chunks} chunks
                approx_file_size = len(records) // num_chunks
                splitted_records = []
                for i in range(num_chunks):
                    start = i * approx_file_size
                    if i != num_chunks -1:
                        end = (i + 1) * approx_file_size
                    else:
                        end = len(records)
                    splitted_records.append(records[start:end])

                # sanity check
                assert len(splitted_records) == num_chunks
                assert sum([len(x) for x in splitted_records]) == len(records)

                # save splitted records to files
                splitted_records = ["".join(record) for record in splitted_records]
                for i, record in enumerate(splitted_records):
                    with open(f"{chunk_dir}/{input_filename}-{str(i).zfill(len(str(num_chunks)))}.fa", "w") as f:
                        f.write(record)

                # merge splitted records to the original file for sanity check
                merged_fa = "".join(splitted_records)
                with open(f"{merged_dir}/{input_filename}.merged.fa", "w") as f:
                    f.write(merged_fa)

                # load merged and original fasta and compare them character by character
                with open(f"{merged_dir}/{input_filename}.merged.fa", "r") as f:
                    merged_file = "".join(f.readlines())
                with open(f"{input_dir}/{input_filename}.fa", "r") as f:
                    original_file = "".join(f.readlines())

                assert original_file == merged_file
