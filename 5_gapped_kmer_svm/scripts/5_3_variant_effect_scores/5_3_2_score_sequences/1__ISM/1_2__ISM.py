import pdb 
import os 
import os.path as osp 

import pandas as pd 
import numpy as np 

from Bio import SeqIO
from tqdm import tqdm


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

# cdSNP gkmpredict score basepath
cdsnp_score_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNP_variant_effect_scores/ISM")
# snATAC cluster names
clusters  = [x for x in os.listdir(cdsnp_score_basepath)]
# CV fold names
fold_ids = os.listdir(cdsnp_score_basepath[0])

# cdSNPs basepath
cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")


# SNP sequence path
seq_filepath_base = f"/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/{analysis_type}"
# gkmpredict target path
gkmpredict_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmpredict"
# ISM score path
ism_filepath_base = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/ISM"


# Step 1: use gkmpredict results to create ISM score folder
seq_len = "201bp"
effect_snp_type = "a1"
non_effect_snp_type = "a2"
seq_types = ["original", "shuffled"]

for cluster in zip(clusters):
    for seq_type in seq_types:
        interim_filepath = osp.join(*[ism_filepath_base, broad_cluster, cell_type, seq_len, seq_type])

        for fold_id in range(1, num_folds+1):

            # load effect and non-effect allele scores
            if not osp.exists(
                    osp.join(*[ 
                        cdsnp_score_basepath,
                        cluster,
                        fold_id,
                        seq_len,
                        seq_type
                    ])
                ):
                print("Skipping {} {} {} {} because it doesn't exist".format(
                    cluster,
                    fold_id,
                    seq_len,
                    seq_type
                ))
                continue
            
            # effect allele
            effect_allele_scores = \
                pd.read_table(
                    osp.join(*[ 
                        cdsnp_score_basepath,
                        cluster,
                        fold_id,
                        seq_len,
                        seq_type,
                        f"cdsnps.unique.{seq_len}.{seq_type}.{effect_snp_type}.gkmpredict"
                    ]),
                    header=None,
                    names=["snp_identifier", "effect_allele_score"]
                )
            # non-effect allele
            non_effect_allele_scores = \
                pd.read_table(
                    osp.join(*[ 
                        cdsnp_score_basepath,
                        cluster,
                        fold_id,
                        seq_len,
                        seq_type,
                        f"cdsnps.unique.{seq_len}.{seq_type}.{non_effect_snp_type}.gkmpredict"
                    ]),
                    header=None,
                    names=["snp_identifier", "effect_allele_score"]
                )

            # create new table for storing allele specific gkmpredict and ISM scores
            ism_scores = effect_allele_scores.copy()
            ism_scores["non_effect_allele_score"] = non_effect_allele_scores["non_effect_allele_score"]

            # calculate ISM scores
            ism_scores["ism_score"] = ism_scores["effect_allele_score"] - ism_scores["non_effect_allele_score"]

            outfile_path = osp.join(*[ 
                        cdsnp_score_basepath,
                        cluster,
                        fold_id,
                        seq_len,
                        seq_type,
                        f"cdsnps.unique.{seq_len}.{seq_type}.{effect_snp_type}_{non_effect_snp_type}.ism"
                    ])
            ism_scores.to_csv(outfile_path, sep="\t", index=False)

print("Script finished!")
