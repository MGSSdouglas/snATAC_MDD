import pdb
import os 
import sys 
import os.path as osp


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

id = "5_3_2_score_sequences/1__ISM"
temp_storage_basepath = osp.join(GKMSVM_TMP_DIR, id)
os.makedirs(temp_storage_basepath, exists_ok=True)


## prepare jobs

# trained gapped-kmer svm model path
models_dir = GKMSVM_MODEL_DIR
# snATAC cluster names
clusters  = [x for x in os.listdir(data_dir)]
# CV fold names
fold_ids = os.listdir(clusters[0])
# gkmsvm model script dir
model_scripts_dir = osp.join(GKMSVM_BIN_DIR, "lsgkm/bin")

# cdSNPs basepath
cdsnp_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNPs")
# cdSNP output score basepath
cdsnp_score_basepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "cdSNP_variant_effect_scores/ISM")
os.makedirs(cdsnp_score_basepath, exists_ok=True)

# osp.join(*[cdsnp_basepath, cluster, f"{seqlen}bp_seqs"])
# f"cdsnps.unique.{seqlen}bp.a1.original.fa"
# f"cdsnps.unique.{seqlen}bp.a2.original.fa"
# f"cdsnps.unique.{seqlen}bp.a1.shuffled.fa"
# f"cdsnps.unique.{seqlen}bp.a2.shuffled.fa"

jobs = []
for cluster in clusters:
    for fold_id in fold_ids:
        if not osp.exists(f"{predict_dir}/{cluster}/{fold_id}"):
            continue
        else:
            jobs.append({
                "Name": f"{cluster}__{fold_id}__ISM",
                "cluster": cluster,
                "fold_id": fold_id,
                "model_path": osp.join(models_dir, cluster, fold_id, f"{cluster}.{fold_id}.model.txt"),
                "model_scripts_dir": model_scripts_dir,
                "cdsnp_dir": cdsnp_basepath,
                "cdsnp_score_dir": cdsnp_score_basepath,
            })

# # # DEBUG
# job_params = [job_params[0]]