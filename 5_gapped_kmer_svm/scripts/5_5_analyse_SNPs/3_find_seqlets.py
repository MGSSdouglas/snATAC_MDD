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
from tqdm import tqdm
from copy import deepcopy
from pprint import pprint
from statsmodels.graphics.gofplots import qqplot_2samples
from statsmodels.distributions.empirical_distribution import ECDF


"""
Configuration
"""

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"
os.makedirs(temp_data_dir, exist_ok=True)

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
output_plot_dir = f"{output_dir}/plots"
os.makedirs(output_data_dir, exist_ok=True)

null_sequences_count = 10

def aggregate_v1(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    if len(np.unique(x)) != 1:
        return ",".join([str(y) for y in x])
    else:
        return [str(y) for y in x][0]

def aggregate_v2(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    return ",".join([str(y) for y in x])

# read latest pipeline_data structure 
with open(f"{temp_data_dir}/merged_data.v5.pkl", "rb") as f:
    pipeline_data = pickle.load(f)

################################################################
# Find 51bp sequences and gkmexplain scores centered at the SNPs
################################################################
middle_pos = 201 // 2
offset = 25

for cell_type in pipeline_data:
    
    
    """
    ['original_snp_tables',
    'shuffled_snp_tables',
    'original_201bp_a1_allele_sequence',
    'original_201bp_a2_allele_sequence',
    'shuffled_201bp_a1_allele_sequence',
    'shuffled_201bp_a2_allele_sequence',
    'original_gkmexplain_ve_scores',
    'shuffled_gkmexplain_ve_scores',
    'original_201bp_gkmexplain_a1_allele_scores',
    'original_201bp_gkmexplain_a2_allele_scores',
    'shuffled_201bp_gkmexplain_a1_allele_scores',
    'shuffled_201bp_gkmexplain_a2_allele_scores',
    'original_deltaSVM_ve_scores',
    'shuffled_deltaSVM_ve_scores',
    'original_ISM_ve_scores',
    'shuffled_ISM_ve_scores',
    'null_ISM_ve_score_related',
    'null_gkmexplain_ve_score_related',
    'null_deltaSVM_ve_score_related']
    """

    # original
    pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"] = \
        pipeline_data[cell_type]["original_201bp_gkmexplain_a1_allele_scores"][:, middle_pos-offset:middle_pos+offset+1, :]

    pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"] = \
        pipeline_data[cell_type]["original_201bp_gkmexplain_a2_allele_scores"][:, middle_pos-offset:middle_pos+offset+1, :]

    # shuffled
    pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"] = \
        pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a1_allele_scores"][:, middle_pos-offset:middle_pos+offset+1, :]

    pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"] = \
        pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a2_allele_scores"][:, middle_pos-offset:middle_pos+offset+1, :]

    # sanity check (1)
    assert pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"].shape[1] == 51
    assert pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"].shape[1] == 51
    assert pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"].shape[1] == 51
    assert pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"].shape[1] == 51

    # """ VERY IMPORTANT NOTE: May need to remove .upper() for the whole script to work """

    # original
    pipeline_data[cell_type]["original_51bp_a1_allele_sequence"] = \
        np.array(
            [list(x)[middle_pos-offset:middle_pos+offset+1] for x in pipeline_data[cell_type]["original_201bp_a1_allele_sequence"]]
        )
    pipeline_data[cell_type]["original_51bp_a2_allele_sequence"] = \
        np.array(
            [list(x)[middle_pos-offset:middle_pos+offset+1] for x in pipeline_data[cell_type]["original_201bp_a2_allele_sequence"]]
        )

    # shuffled
    pipeline_data[cell_type]["shuffled_51bp_a1_allele_sequence"] = \
        np.array(
            [list(x)[middle_pos-offset:middle_pos+offset+1] for x in pipeline_data[cell_type]["shuffled_201bp_a1_allele_sequence"]]
        )
    pipeline_data[cell_type]["shuffled_51bp_a2_allele_sequence"] = \
        np.array(
            [list(x)[middle_pos-offset:middle_pos+offset+1] for x in pipeline_data[cell_type]["shuffled_201bp_a2_allele_sequence"]]
        )
    
    
    # sanity check (2)
    assert pipeline_data[cell_type]["original_51bp_a1_allele_sequence"].shape[1] == 51
    assert pipeline_data[cell_type]["original_51bp_a2_allele_sequence"].shape[1] == 51
    assert pipeline_data[cell_type]["shuffled_51bp_a1_allele_sequence"].shape[1] == 51
    assert pipeline_data[cell_type]["shuffled_51bp_a2_allele_sequence"].shape[1] == 51

    ### sanity check (3)
    dim2base = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

    ## original_201bp_a1_allele_sequence vs original_201bp_gkmexplain_a1_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["original_201bp_gkmexplain_a1_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare 
    assert np.all(tempp == pipeline_data[cell_type]["original_201bp_a1_allele_sequence"])

    ## original_51bp_a1_allele_sequence vs original_51bp_gkmexplain_a1_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["original_51bp_a1_allele_sequence"])

    ## original_201bp_a2_allele_sequence vs original_201bp_gkmexplain_a2_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["original_201bp_gkmexplain_a2_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["original_201bp_gkmexplain_a2_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["original_201bp_gkmexplain_a2_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["original_201bp_a2_allele_sequence"])

    ## original_51bp_a2_allele_sequence vs original_51bp_gkmexplain_a2_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["original_51bp_a2_allele_sequence"])

    ## shuffled_201bp_a1_allele_sequence vs shuffled_201bp_gkmexplain_a1_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a1_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a1_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a1_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["shuffled_201bp_a1_allele_sequence"])

    ## shuffled_51bp_a1_allele_sequence vs shuffled_51bp_gkmexplain_a1_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["shuffled_51bp_a1_allele_sequence"])

    ## shuffled_201bp_a2_allele_sequence vs shuffled_201bp_gkmexplain_a2_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a2_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a2_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a2_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["shuffled_201bp_a2_allele_sequence"])

    ## shuffled_51bp_a2_allele_sequence vs shuffled_51bp_gkmexplain_a2_allele_scores

    # get index of absolute max of second dimension of pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"]
    temp = np.argmax(np.abs(pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"]), axis=-1).tolist()
    # get base of absolute max of second dimension of pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"]
    tempp = np.array([[dim2base[x] for x in row] for row in temp])
    # compare
    assert np.all(tempp == pipeline_data[cell_type]["shuffled_51bp_a2_allele_sequence"])


print("All tests seq vs seq based on gkmexplain scores are passed!")
    

"""
Current keys of pipeline data:

['original_snp_tables',
 'shuffled_snp_tables',
 'original_201bp_a1_allele_sequence',
 'original_201bp_a2_allele_sequence',
 'shuffled_201bp_a1_allele_sequence',
 'shuffled_201bp_a2_allele_sequence',
 'original_gkmexplain_ve_scores',
 'shuffled_gkmexplain_ve_scores',
 'original_201bp_gkmexplain_a1_allele_scores',
 'original_201bp_gkmexplain_a2_allele_scores',
 'shuffled_201bp_gkmexplain_a1_allele_scores',
 'shuffled_201bp_gkmexplain_a2_allele_scores',
 'original_deltaSVM_ve_scores',
 'shuffled_deltaSVM_ve_scores',
 'original_ISM_ve_scores',
 'shuffled_ISM_ve_scores',
 'null_ISM_ve_score_related',
 'null_gkmexplain_ve_score_related',
 'null_deltaSVM_ve_score_related',

    The following keys are newly added

 'original_51bp_gkmexplain_a1_allele_scores',
 'original_51bp_gkmexplain_a2_allele_scores',
 'shuffled_51bp_gkmexplain_a1_allele_scores',
 'shuffled_51bp_gkmexplain_a2_allele_scores',
 'original_51bp_a1_allele_sequence',
 'original_51bp_a2_allele_sequence',
 'shuffled_51bp_a1_allele_sequence',
 'shuffled_51bp_a2_allele_sequence']
"""

#############################################################
## Find active alleles for (SNP, cell_type) pairs
#############################################################
for cell_type in pipeline_data:

    """
    Active allele is defined as the allele having higher sum of non-negative
    importance scores in the 51bp region centered at the SNP 
    """
    # original
    # a1
    scores = pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]
    a1_non_negative_pos_sum = np.where(scores>=0,scores,0).sum((-1,-2))
    # a2
    scores = pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"]
    a2_non_negative_pos_sum = np.where(scores>=0,scores,0).sum((-1,-2))
    # if a1 sum is greater than or equal to a2 sum, then a1 is active
    # if a2 sum is greater than a1 sum, then a2 is active
    active_alleles = np.where(a1_non_negative_pos_sum>=a2_non_negative_pos_sum, 'A1', 'A2')
    # save active alleles
    pipeline_data[cell_type]["original_snp_tables"]["active_allele"] = active_alleles
    # create active allele score array 
    active_allele_score_array = []
    active_allele_seq_array = []
    for ix, active_allele in enumerate(active_alleles):
        if active_allele == "A1":
            score_source = "original_51bp_gkmexplain_a1_allele_scores"
            seq_source = "original_51bp_a1_allele_sequence"
        elif active_allele == "A2":
            score_source = "original_51bp_gkmexplain_a2_allele_scores"
            seq_source = "original_51bp_a2_allele_sequence"
        # store active allele sequence and scores
        active_allele_score_array.append(
            pipeline_data[cell_type][score_source][ix]
        )
        active_allele_seq_array.append(
            pipeline_data[cell_type][seq_source][ix]
        )
    active_allele_score_array = np.array(active_allele_score_array)
    active_allele_seq_array = np.array(active_allele_seq_array)
    # save active allele score array
    pipeline_data[cell_type]["original_51bp_gkmexplain_active_allele_scores"] = active_allele_score_array
    pipeline_data[cell_type]["original_51bp_active_allele_sequence"] = active_allele_seq_array

    # shuffled
    # a1
    scores = pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"]
    a1_non_negative_pos_sum = np.where(scores>=0,scores,0).sum((-1,-2))
    # a2
    scores = pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"]
    a2_non_negative_pos_sum = np.where(scores>=0,scores,0).sum((-1,-2))
    # if a1 sum is greater than or equal to a2 sum, then a1 is active
    # if a2 sum is greater than a1 sum, then a2 is active
    active_alleles = np.where(a1_non_negative_pos_sum>=a2_non_negative_pos_sum, 'A1', 'A2')
    # save active alleles
    pipeline_data[cell_type]["shuffled_snp_tables"]["active_allele"] = active_alleles
    # create active allele score array 
    active_allele_score_array = []
    active_allele_seq_array = []
    for ix, active_allele in enumerate(active_alleles):
        if active_allele == "A1":
            score_source = "shuffled_51bp_gkmexplain_a1_allele_scores"
            seq_source = "shuffled_51bp_a1_allele_sequence"
        elif active_allele == "A2":
            score_source = "shuffled_51bp_gkmexplain_a2_allele_scores"
            seq_source = "shuffled_51bp_a2_allele_sequence"
        # store active allele sequence and scores
        active_allele_score_array.append(
            pipeline_data[cell_type][score_source][ix]
        )
        active_allele_seq_array.append(
            pipeline_data[cell_type][seq_source][ix]
        )
    active_allele_score_array = np.array(active_allele_score_array)
    active_allele_seq_array = np.array(active_allele_seq_array)
    # save active allele score array
    pipeline_data[cell_type]["shuffled_51bp_gkmexplain_active_allele_scores"] = active_allele_score_array
    pipeline_data[cell_type]["shuffled_51bp_active_allele_sequence"] = active_allele_seq_array

# annotate shuffled snp table with statistical significance
for cell_type in pipeline_data:
    # extract (rsid, statistical significance) pairs from original snp tables
    original_snp_table = pipeline_data[cell_type]["original_snp_tables"]
    rsids = original_snp_table["Name"].tolist()
    statistical_significance = original_snp_table["ssSNP"].tolist()  
    rsid2ss = dict(zip(rsids, statistical_significance))
    # annotate snps in shuffled table with statistical significance
    shuffled_snp_table = pipeline_data[cell_type]["shuffled_snp_tables"]
    shuffled_snp_table["ssSNP"] = shuffled_snp_table["Name"].map(rsid2ss)
    pipeline_data[cell_type]["shuffled_snp_tables"] = shuffled_snp_table

######################################################################################
## Find 97.5% percentile of single bp gkmexplain scores on 51bp shuffled data
## NOTE: this step is done for all candidate snps instead of statistically significant snps
######################################################################################
percentile = 97.5
for cell_type in pipeline_data:
    # create tf motif related data structure 
    pipeline_data[cell_type]["pipeline_tf_motif_related"] = {}

    # extract shuffled 201bp sequence single bp gkmexplain scores 
    shuffled_scores = []
    # a1
    scores = pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"]
    shuffled_scores += scores[scores>=0].flatten().tolist()
    # a2
    scores = pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"]
    shuffled_scores += scores[scores>=0].flatten().tolist()
    # save
    pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_51bp_gkmexplain_scores"] = shuffled_scores

    # find 97.5% percentile of single bo gkmexplain scores on 201bp shuffled data
    percentile_score = np.percentile(shuffled_scores, percentile, method='lower')
    pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_51bp_gkmexplain_97.5_thr"] = percentile_score

#################################
## Find seqlet for active alleles
#################################

middle_pos_51bp = 25

dim2base = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
}

def find_seqlet(active_allele_scores,
                percentile_threshold, 
                snp_pos,
                min_seqlet_length=7):
    """This function determines active allele seqlet
    """

    # initialize
    seqlet = f"{dim2base[np.argmax(np.absolute(active_allele_scores[snp_pos,:]))]}"
    seqlet_start_pos = snp_pos
    seqlet_end_pos = snp_pos
    seqlet_length = 1

    # iterate
    tolerance = 2
    start_pos_tolerance = tolerance
    end_pos_tolerance = tolerance
    current_dir = "upstream"
    while start_pos_tolerance > 0 or end_pos_tolerance > 0:

        # reached to bounds of 51bp sequence -> no extend
        if seqlet_start_pos == 0 or seqlet_end_pos == 50:
            break

        # extend
        elif current_dir == "upstream":
            if start_pos_tolerance > 0:
                # extend upstream if condition is satisfied
                proposed_start_pos = seqlet_start_pos - 1
                cur_score = active_allele_scores[proposed_start_pos,:].sum()
                # current score greater than null background threshold
                if cur_score >= percentile_threshold:
                    seqlet_start_pos = proposed_start_pos
                    seqlet_length += 1
                    start_pos_tolerance = tolerance
                    # update seqlet
                    seqlet = dim2base[np.argmax(np.absolute(active_allele_scores[seqlet_start_pos,:]))] + seqlet
                # current score less than null background threshold
                else:
                    seqlet_start_pos = proposed_start_pos
                    seqlet_length += 1
                    # still update seqlet but decrease tolerance
                    seqlet = dim2base[np.argmax(np.absolute(active_allele_scores[seqlet_start_pos,:]))] + seqlet
                    start_pos_tolerance -= 1
            # invert direction
            current_dir = "downstream"

        elif current_dir == "downstream":
            if end_pos_tolerance > 0:
                # extend downstream if condition is satisfied
                proposed_end_pos = seqlet_end_pos + 1
                cur_score = active_allele_scores[proposed_end_pos,:].sum()
                # current score greater than null background threshold
                if cur_score >= percentile_threshold:
                    seqlet_end_pos = proposed_end_pos
                    seqlet_length += 1
                    end_pos_tolerance = tolerance
                    # update seqlet
                    seqlet = seqlet + dim2base[np.argmax(np.absolute(active_allele_scores[seqlet_end_pos,:]))]
                # current score less than null background threshold
                else:
                    seqlet_end_pos = proposed_end_pos
                    seqlet_length += 1
                    # still update seqlet but decrease tolerance
                    seqlet = seqlet + dim2base[np.argmax(np.absolute(active_allele_scores[seqlet_end_pos,:]))]
                    end_pos_tolerance -= 1
            # invert direction
            current_dir = "upstream"

    if len(seqlet) != seqlet_length:
        pdb.set_trace()
        print("Should not be here (length,1)")

    # adjust start and end 
    seqlet_start_pos = seqlet_start_pos + (tolerance - start_pos_tolerance)
    seqlet_end_pos = seqlet_end_pos - (tolerance - end_pos_tolerance)
    seqlet_length = seqlet_end_pos - seqlet_start_pos + 1
    if tolerance - end_pos_tolerance == 0:
        seqlet = seqlet[(tolerance - start_pos_tolerance):]
    else:
        seqlet = seqlet[(tolerance - start_pos_tolerance):-(tolerance - end_pos_tolerance)]

    if len(seqlet) != seqlet_length:
        pdb.set_trace()
        print("Should not be here (length,2)")

    # if seqlet is too short
    # alternatingly extend starting from current direction
    # stop when seqlet is long enough
    if seqlet_length < min_seqlet_length:
        #  until minimum length is reached
        loop_ctr = 0
        while seqlet_length < min_seqlet_length:
            if current_dir == "upstream" and seqlet_start_pos > 0:
                seqlet_start_pos -= 1
                cur_score = active_allele_scores[seqlet_start_pos,:].sum()
                seqlet = dim2base[np.argmax(np.absolute(active_allele_scores[seqlet_start_pos,:]))] + seqlet
                current_dir = "downstream"
            elif current_dir == "downstream" and seqlet_end_pos < 50:
                seqlet_end_pos += 1
                cur_score = active_allele_scores[seqlet_end_pos,:].sum()
                seqlet = seqlet + dim2base[np.argmax(np.absolute(active_allele_scores[seqlet_end_pos,:]))]
                current_dir = "upstream"
            loop_ctr += 1
            seqlet_length += 1
            if loop_ctr > min_seqlet_length:
                pdb.set_trace()
                print("Should not be here (loop)")

    assert seqlet_length == len(seqlet) >= min_seqlet_length

    return (seqlet, seqlet_start_pos, seqlet_end_pos)


for cell_type in pipeline_data:
    # get 97.5 percentile threshold
    percentile_threshold = \
        pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_51bp_gkmexplain_97.5_thr"]
    
    # original
    # get active allele 51bp gkmexplain scores
    active_allele_scores = \
        pipeline_data[cell_type]["original_51bp_gkmexplain_active_allele_scores"]

    # convert first dimension of the scores to a python list 
    active_allele_scores = \
        [active_allele_scores[i] for i in range(active_allele_scores.shape[0])]
    # apply `find_seqlet` function to each element in the list
    seqlets, seqlet_starts, seqlet_ends = [], [], []
    for scores in active_allele_scores:
        # extract seqlet
        seqlet, seqlet_start, seqlet_end = \
            find_seqlet(scores, percentile_threshold, middle_pos_51bp)
        # sanity check (3): seqlets should have at least 7bp length
        assert len(seqlet) >= 7 and seqlet_end - seqlet_start + 1 >= 7, "Sanity check (3) is failed"
        # record seqlet details
        seqlets.append(seqlet)
        seqlet_starts.append(seqlet_start)
        seqlet_ends.append(seqlet_end)
    # add seqlet details to the pipeline 
    pipeline_data[cell_type]["original_snp_tables"]["seqlet"] = seqlets
    pipeline_data[cell_type]["original_snp_tables"]["seqlet_start"] = seqlet_starts
    pipeline_data[cell_type]["original_snp_tables"]["seqlet_end"] = seqlet_ends

    # shuffled
    # get active allele 51bp gkmexplain scores
    active_allele_scores = \
        pipeline_data[cell_type]["shuffled_51bp_gkmexplain_active_allele_scores"]
    # convert first dimension of the scores to a python list 
    active_allele_scores = \
        [active_allele_scores[i] for i in range(active_allele_scores.shape[0])]
    # apply `find_seqlet` function to each element in the list
    seqlets, seqlet_starts, seqlet_ends = [], [], []
    for scores in active_allele_scores:
        # extract seqlet
        seqlet, seqlet_start, seqlet_end = \
            find_seqlet(scores, percentile_threshold, middle_pos_51bp)
        # sanity check (4): seqlets should have at least 7bp length
        assert len(seqlet) >= 7 and seqlet_end - seqlet_start + 1 >= 7, "Sanity check (4) is failed"
        # record seqlet details
        seqlets.append(seqlet)
        seqlet_starts.append(seqlet_start)
        seqlet_ends.append(seqlet_end)
    # add seqlet details to the pipeline 
    pipeline_data[cell_type]["shuffled_snp_tables"]["seqlet"] = seqlets
    pipeline_data[cell_type]["shuffled_snp_tables"]["seqlet_start"] = seqlet_starts
    pipeline_data[cell_type]["shuffled_snp_tables"]["seqlet_end"] = seqlet_ends

### sanity check:

for cell_type in pipeline_data:

    ## original

    snp_table = pipeline_data[cell_type]["original_snp_tables"]

    # seqlet details 
    seqlet_start_pos = snp_table["seqlet_start"].astype(int).tolist()
    seqlet_end_pos = snp_table["seqlet_end"].astype(int).tolist()
    seqlet = snp_table["seqlet"].tolist()

    # get active allele 51bp gkmexplain scores
    active_allele_scores = pipeline_data[cell_type]["original_51bp_gkmexplain_active_allele_scores"].tolist()
    active_allele_seqlet_scores = [active_allele_scores[i][seqlet_start_pos[i]:seqlet_end_pos[i]+1] for i in range(len(active_allele_scores))]

    # get active allele seqlets
    active_allele_seqs = pipeline_data[cell_type]["original_51bp_active_allele_sequence"].tolist()
    active_allele_seqlets = ["".join(active_allele_seqs[i]).upper()[seqlet_start_pos[i]:seqlet_end_pos[i]+1] for i in range(len(active_allele_seqs))]

    # active_allele_snp_pos = [x[middle_pos_51bp] for x in active_allele_seqs]

    try:
        assert [True if x == y else False for x, y in zip(active_allele_seqlets, seqlet)]
    except:
        print(f"{cell_type}\t seqlets != active_allele_seqlets")
        pdb.set_trace()

##################################################################
## 1 - Calculate `seqlet score` and `seqlet SNR` for A1 and A2 alleles 
## 2 - Calculate `prominence` and `magnitude` scores from `seqlet score` and `seqlet SNR` of each allele
##################################################################
for cell_type in pipeline_data:

    ### original
    original_snp_table = pipeline_data[cell_type]["original_snp_tables"]

    # a1
    original_a1_allele_51bp_gkmexplain_scores = \
        pipeline_data[cell_type]["original_51bp_gkmexplain_a1_allele_scores"]
    original_a1_allele_201bp_gkmexplain_scores = \
        pipeline_data[cell_type]["original_201bp_gkmexplain_a1_allele_scores"]
    # seqlet score and seqlet signal-to-noise ratio
    a1_allele_seqlet_scores, a1_allele_seqlet_snrs = [], []
    for ix in range(original_a1_allele_201bp_gkmexplain_scores.shape[0]):
        # prepare required variables
        gkmexplain_scores_51bp = original_a1_allele_51bp_gkmexplain_scores[ix,:]
        gkmexplain_scores_201bp = original_a1_allele_201bp_gkmexplain_scores[ix,:]
        active_allele_seqlet_start = original_snp_table["seqlet_start"].iloc[ix]
        active_allele_seqlet_end = original_snp_table["seqlet_end"].iloc[ix]
        # calculate seqlet score: sum non-negative entries of 51bp score for positions overlapping the seqlet
        seqlet_score = gkmexplain_scores_51bp[seqlet_start:seqlet_end+1]
        seqlet_score = np.sum(seqlet_score[seqlet_score>0]) 
        # calculate seqlet SNR: seqlet score / sum non-negative entries of 201bp score
        denom = np.sum(gkmexplain_scores_201bp[gkmexplain_scores_201bp>0])
        if denom == 0: # to avoid divide-by-zero error
            seqlet_snr = np.inf
        else:
            seqlet_snr = seqlet_score / denom
        # store these scores
        a1_allele_seqlet_scores.append(seqlet_score)
        a1_allele_seqlet_snrs.append(seqlet_snr)
    # record these values to original snp table 
    original_snp_table["a1_allele_seqlet_score"] = a1_allele_seqlet_scores
    original_snp_table["a1_allele_seqlet_snr"] = a1_allele_seqlet_snrs

    # a2
    original_a2_allele_51bp_gkmexplain_scores = \
        pipeline_data[cell_type]["original_51bp_gkmexplain_a2_allele_scores"]
    original_a2_allele_201bp_gkmexplain_scores = \
        pipeline_data[cell_type]["original_201bp_gkmexplain_a2_allele_scores"]
    # seqlet score and seqler signal-to-noise ratio
    a2_allele_seqlet_scores, a2_allele_seqlet_snrs = [], []
    for ix in range(original_a2_allele_201bp_gkmexplain_scores.shape[0]):
        # prepare required variables
        gkmexplain_scores_51bp = original_a2_allele_51bp_gkmexplain_scores[ix,:]
        gkmexplain_scores_201bp = original_a2_allele_201bp_gkmexplain_scores[ix,:]
        active_allele_seqlet_start = original_snp_table["seqlet_start"].iloc[ix]
        active_allele_seqlet_end = original_snp_table["seqlet_end"].iloc[ix]
        # calculate seqlet score: sum non-negative entries of 51bp score for positions overlapping the seqlet
        seqlet_score = gkmexplain_scores_51bp[seqlet_start:seqlet_end+1]
        seqlet_score = np.sum(seqlet_score[seqlet_score>0]) 
        # calculate seqlet SNR: seqlet score / sum non-negative entries of 201bp score
        denom = np.sum(gkmexplain_scores_201bp[gkmexplain_scores_201bp>0])
        if denom == 0:
            seqlet_snr = np.inf
        else:
            seqlet_snr = seqlet_score / denom
        # store these scores
        a2_allele_seqlet_scores.append(seqlet_score)
        a2_allele_seqlet_snrs.append(seqlet_snr)
    # record these values to original snp table
    original_snp_table["a2_allele_seqlet_score"] = a2_allele_seqlet_scores
    original_snp_table["a2_allele_seqlet_snr"] = a2_allele_seqlet_snrs

    # calculate magnitude score: a1_seqlet_score - a2_seqlet_score
    original_snp_table["magnitude_score"] = \
        original_snp_table["a1_allele_seqlet_score"] - original_snp_table["a2_allele_seqlet_score"]
    # calculate prominence score: a1_seqlet_snr - a2_seqlet_snr
    original_snp_table["prominence_score"] = \
        original_snp_table["a1_allele_seqlet_snr"] - original_snp_table["a2_allele_seqlet_snr"]

    # make changes to pipeline data structure
    pipeline_data[cell_type]["original_snp_tables"] = original_snp_table

    ### shuffled 
    shuffled_snp_table = pipeline_data[cell_type]["shuffled_snp_tables"]

    # a1
    shuffled_a1_allele_51bp_gkmexplain_scores = \
        pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a1_allele_scores"]
    shuffled_a1_allele_201bp_gkmexplain_scores = \
        pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a1_allele_scores"]
    # seqlet score and seqler signal-to-noise ratio
    a1_allele_seqlet_scores, a1_allele_seqlet_snrs = [], []
    for ix in range(shuffled_a1_allele_201bp_gkmexplain_scores.shape[0]):
        # prepare required variables
        gkmexplain_scores_51bp = shuffled_a1_allele_51bp_gkmexplain_scores[ix,:]
        gkmexplain_scores_201bp = shuffled_a1_allele_201bp_gkmexplain_scores[ix,:]
        active_allele_seqlet_start = shuffled_snp_table["seqlet_start"].iloc[ix]
        active_allele_seqlet_end = shuffled_snp_table["seqlet_end"].iloc[ix]
        # calculate seqlet score: sum non-negative entries of 51bp score for positions overlapping the seqlet
        seqlet_score = gkmexplain_scores_51bp[seqlet_start:seqlet_end+1]
        seqlet_score = np.sum(seqlet_score[seqlet_score>0]) 
        # calculate seqlet SNR: seqlet score / sum non-negative entries of 201bp score
        denom = np.sum(gkmexplain_scores_201bp[gkmexplain_scores_201bp>0])
        if denom == 0:
            seqlet_snr = np.inf
        else:
            seqlet_snr = seqlet_score / denom
        # store these scores
        a1_allele_seqlet_scores.append(seqlet_score)
        a1_allele_seqlet_snrs.append(seqlet_snr)
    # record these values to shuffled snp table
    shuffled_snp_table["a1_allele_seqlet_score"] = a1_allele_seqlet_scores
    shuffled_snp_table["a1_allele_seqlet_snr"] = a1_allele_seqlet_snrs

    # a2
    shuffled_a2_allele_51bp_gkmexplain_scores = \
        pipeline_data[cell_type]["shuffled_51bp_gkmexplain_a2_allele_scores"]
    shuffled_a2_allele_201bp_gkmexplain_scores = \
        pipeline_data[cell_type]["shuffled_201bp_gkmexplain_a2_allele_scores"]
    # seqlet score and seqler signal-to-noise ratio
    a2_allele_seqlet_scores, a2_allele_seqlet_snrs = [], []
    for ix in range(shuffled_a2_allele_201bp_gkmexplain_scores.shape[0]):
        # prepare required variables
        gkmexplain_scores_51bp = shuffled_a2_allele_51bp_gkmexplain_scores[ix,:]
        gkmexplain_scores_201bp = shuffled_a2_allele_201bp_gkmexplain_scores[ix,:]
        active_allele_seqlet_start = shuffled_snp_table["seqlet_start"].iloc[ix]
        active_allele_seqlet_end = shuffled_snp_table["seqlet_end"].iloc[ix]
        # calculate seqlet score: sum non-negative entries of 51bp score for positions overlapping the seqlet
        seqlet_score = gkmexplain_scores_51bp[seqlet_start:seqlet_end+1]
        seqlet_score = np.sum(seqlet_score[seqlet_score>0]) 
        # calculate seqlet SNR: seqlet score / sum non-negative entries of 201bp score
        denom = np.sum(gkmexplain_scores_201bp[gkmexplain_scores_201bp>0])
        if denom == 0:
            seqlet_snr = np.inf
        else:
            seqlet_snr = seqlet_score / denom
        # store these scores
        a2_allele_seqlet_scores.append(seqlet_score)
        a2_allele_seqlet_snrs.append(seqlet_snr)
    # record these values to shuffled snp table
    shuffled_snp_table["a2_allele_seqlet_score"] = a2_allele_seqlet_scores
    shuffled_snp_table["a2_allele_seqlet_snr"] = a2_allele_seqlet_snrs

    # calculate magnitude score: a1_seqlet_score - a2_seqlet_score
    shuffled_snp_table["magnitude_score"] = \
        shuffled_snp_table["a1_allele_seqlet_score"] - shuffled_snp_table["a2_allele_seqlet_score"]
    # calculate prominence score: a1_seqlet_snr - a2_seqlet_snr
    shuffled_snp_table["prominence_score"] = \
        shuffled_snp_table["a1_allele_seqlet_snr"] - shuffled_snp_table["a2_allele_seqlet_snr"]

    # make changes to pipeline data structure
    pipeline_data[cell_type]["shuffled_snp_table"] = shuffled_snp_table

############################################################
## Find null distribution of prominence and magnitude scores
############################################################

for cell_type in pipeline_data:
    # get shuffled snp table
    shuffled_snp_table = pipeline_data[cell_type]["shuffled_snp_table"]
    # get shuffled magnitude and prominence scores
    shuffled_magnitude_scores = shuffled_snp_table["magnitude_score"].tolist()
    shuffled_prominence_scores = shuffled_snp_table["prominence_score"].tolist()
    # add negative of each score to the null distribution to make the distribution symmetric
    # this is to make the null distribution symmetric around 0
    shuffled_magnitude_scores += [-x for x in shuffled_magnitude_scores]
    shuffled_prominence_scores += [-x for x in shuffled_prominence_scores]
    
    # compute ECDF of the shuffled scores 
    shufled_magnitude_score_ECDF = ECDF(shuffled_magnitude_scores)
    shufled_prominence_score_ECDF = ECDF(shuffled_prominence_scores)
    # sort magnitude and prominence scores for plotting purposes
    shuffled_magnitude_scores.sort()
    shuffled_prominence_scores.sort()
    # get ECDF values for the sorted scores 
    shuffled_magnitude_score_ECDF_values = shufled_magnitude_score_ECDF(shuffled_magnitude_scores)
    shuffled_prominence_score_ECDF_values = shufled_prominence_score_ECDF(shuffled_prominence_scores)

    # plot magnitude ECDF
    plt.figure()
    plt.plot(shuffled_magnitude_scores, shuffled_magnitude_score_ECDF_values, "b-")
    plt.xlabel(r"Magnitude scores ($x_{obs}$)")
    plt.ylabel(r"Probability $P(X \leq x_{obs})$")
    plt.title("Symmetric null magnitude score ECDF for " + cell_type)
    plt.tight_layout()
    plt.savefig(f"{output_plot_dir}/{cell_type}/null_magnitude_ECDF.png")
    plt.close()

    # plot prominence ECDF
    plt.figure()
    plt.plot(shuffled_prominence_scores, shuffled_prominence_score_ECDF_values, "b-")
    plt.xlabel(r"Prominence scores ($x_{obs}$)")
    plt.ylabel(r"Probability $P(X \leq x_{obs})$")
    plt.title("Symmetric null magnitude score ECDF for " + cell_type)
    plt.tight_layout()
    plt.savefig(f"{output_plot_dir}/{cell_type}/null_prominence_ECDF.png")
    plt.close()

    # save these null distribution detail to pipeline_data
    pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_magnitude_scores"] = \
        shuffled_magnitude_scores
    pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_prominence_scores"] = \
        shuffled_prominence_scores
    pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_magnitude_score_ECDF"] = \
        shufled_magnitude_score_ECDF
    pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_prominence_score_ECDF"] = \
        shufled_prominence_score_ECDF

##############################################
## Calculate p-values for candidate and ssSNPs
##############################################

calculate_p_value = lambda score, ecdf: \
    2 * ecdf(score) if ecdf(score) <= 0.5 else 2 * (1 - ecdf(score))

for cell_type in pipeline_data:
    # get relevant data
    original_snp_table = pipeline_data[cell_type]["original_snp_tables"]

    magnitude_scores = original_snp_table["magnitude_score"].tolist()
    magnitude_score_ecdf = pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_magnitude_score_ECDF"]
    prominence_scores = original_snp_table["prominence_score"].tolist()
    prominence_score_ecdf = pipeline_data[cell_type]["pipeline_tf_motif_related"]["shuffled_prominence_score_ECDF"]
    
    # calculate p-values for magnitude score 
    magnitude_score_pval = []
    for score in magnitude_scores:
        magnitude_score_pval.append(calculate_p_value(score, magnitude_score_ecdf))

    # calculate p-values for prominence score
    prominence_score_pval = []
    for score in prominence_scores:
        prominence_score_pval.append(calculate_p_value(score, prominence_score_ecdf))

    # assign confidence levels to snps using p-values calculated above 
    conf_levels = []
    for m_pval, p_pval in zip(magnitude_score_pval, prominence_score_pval):
        if p_pval < 0.05:
            conf_levels.append("High")
        elif m_pval < 0.05 or p_pval < 0.1:
            conf_levels.append("Medium")
        else:
            conf_levels.append("Low")
    
    # store these values
    original_snp_table["magnitude_pval"] = magnitude_score_pval
    original_snp_table["prominence_pval"] = prominence_score_pval
    original_snp_table["tf_motif_distruption_confidence"] = conf_levels

    # save this table to pipeline_data
    pipeline_data[cell_type]["original_snp_tables"] = original_snp_table

# save the resulting pipeline_data structure 
with open(f"{temp_data_dir}/merged_data.v6.pkl", "wb") as f:
    pickle.dump(pipeline_data, f)

# flatten snp_tables across cell types for eQTL analysis
tables = []
for cell_type in pipeline_data:
    tables.append(pipeline_data[cell_type]["original_snp_tables"])
snp_table = pd.concat(tables, axis=0)
snp_table.to_csv(f"{temp_data_dir}/merged_data.v6.ct_independent.csv", index=False)
    
print("Script finished")
