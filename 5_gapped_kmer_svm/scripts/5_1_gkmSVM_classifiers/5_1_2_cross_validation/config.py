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

id = "5_1_gkmSVM_classifiers/5_1_2_cross_validation"
temp_storage_basepath = osp.join(GKMSVM_TMP_DIR, "5_1_gkmSVM_classifiers/5_1_2_cross_validation")
os.makedirs(temp_storage_basepath, exists_ok=True)

## prepare jobs

# source filepaths
id = "5_1_gkmSVM_classifiers/5_1_1_preprocessing"
data_dir = osp.join(GKMSVM_TMP_DIR, id)
clusters  = [x for x in os.listdir(data_dir)]
fold_ids = os.listdir(clusters[0])

# model script filepaths 
model_scripts_dir = os.environ.get("GKMSVM")

# model storage filepaths
models_dir = GKMSVM_MODEL_DIR
os.makedirs(models_dir, exist_ok=True)

# prediction filepaths
id = "5_1_gkmSVM_classifiers/5_1_2_cross_validation"
predict_dir = osp.join(GKMSVM_TMP_DIR, id)
os.makedirs(predict_dir, exist_ok=True)

# job params
jobs = []
for cluster in clusters:
    for fold_id in fold_ids:
        if osp.exists(f"{predict_dir}/{cluster}/{fold_id}"):
            continue
        else:
            jobs.append({"Name": f"{cluster}_{fold_id}",
                    "cluster": cluster,
                    "fold_id": fold_id,
                    "data_dir": data_dir,
                    "predict_dir": predict_dir,
                    "model_scripts_dir": model_scripts_dir,
                    "models_dir": models_dir})
