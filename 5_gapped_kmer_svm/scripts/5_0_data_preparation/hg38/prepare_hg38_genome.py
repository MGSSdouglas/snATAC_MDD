import os
import re
import pdb
import random
import pickle

import os.path as osp
import numpy as np

from Bio import SeqIO


def check_available(genome_dir, chr_filenames):
    """ 
    This function checks if the files for the specified genome is
        ready for preprocessing
    :@param genome_path: File path to genome
    :@param chr_filenames: A list of file names for the chromosome that 
        should have been downloaded
    """
    assert osp.exists(genome_dir),\
        f"Please run bash script for downloading the requested genome first."

    for filename in chr_filenames:
        assert osp.exists(osp.join(genome_dir, filename)),\
            f"Please run bash script for downloading the requested genome first."


def load_genome_to_dict(genome_dir, chr_filenames):
    """Load the specified genome to a dictionary where keys are chromosome
    names and values are genomic sequence for that chromosome
    :@param genome_path: File path to genome
    :@param chr_filenames: A list of file names for the chromosome that 
        should have been downloaded
    
    NOTE: This function is for reading chromosomes from their respective fasta files

    """

    check_available(genome_dir, chr_filenames)

    chr_filenames = [osp.join(genome_dir, chrom)\
        for chrom in chr_filenames]

    # load fasta files per chromosome
    sequences = {}
    for filename in chr_filenames:
        file = open(filename)
        fasta = SeqIO.parse(file, "fasta")
        for record in fasta:
            name, sequence = record.id, str(record.seq)
            sequences[name] = sequence
        file.close()
    
    return sequences


def load_genome_from_single_fasta_to_dict(fasta_path, verbosity=0):
    '''Load genome from a single fasta file to a dictionary where keys are chromosome

    :param fasta_path: File path to fasta file
    :param verbosity: Verbosity level (if 0 no print, if greater than 0 then print)
    :return sequences: Dictionary of chromosome names corresponding genomic sequences

    NOTE: regular expression pattern in the function is hardcoded to recognize chr[1-22] + chr[X-Y]
    '''

    file = open(fasta_path)
    fasta_sequences = SeqIO.parse(file, "fasta")

    # load chromosomes in a fasta file one by one
    pattern = re.compile(r'^chr[0-9]*[XY]*$')
    sequences = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if pattern.match(name) != None:
            sequences[name] = sequence
            if verbosity > 0:
                print(name, f"sequence length: {len(sequence)}")
                print()
        else:
            print(f"{name} does not match the regex pattern")
            print()

    return sequences

if __name__ == "__main__":

    # get environment variables
    GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
    GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
    GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
    GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
    GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")

    data_id = "hg38"
    data_dir = f"{GKMSVM_RAW_DATA_DIR}/{data_id}"

    fasta_path = f"{data_dir}/hg38.p13.fa"
    dct = load_genome_from_single_fasta_to_dict(fasta_path, verbosity=1)

    save_fname = osp.join(data_dir, "grch38_p13.dict.pkl")
    with open(save_fname, "wb") as f:
        pickle.dump(dct, f)

    print("Script finished")
