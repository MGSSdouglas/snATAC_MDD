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

## prepare jobs

# trained gapped-kmer svm model path
models_dir = GKMSVM_MODEL_DIR
# snATAC cluster names
clusters  = [x for x in os.listdir(data_dir)]
# CV fold names
fold_ids = os.listdir(clusters[0])
# gkmsvm model script dir
model_scripts_dir = osp.join(GKMSVM_BIN_DIR, "lsgkm/bin")


'''
I/O file paths
'''

# SNP sequence path
seq_base_path = f"/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/{analysis_type}"

# deltaSVM target path
deltasvm_score_base_path = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/deltaSVM"

# gkmpredict target path
gkmpredict_score_base_path = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmpredict"

# kmer fasta path 
kmer_fasta_path = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/all_11mers.fa"


'''
Create job parameters
'''

job_cell_types = []
job_broad_clusters = []
job_fold_ids = []
job_model_base_paths = []
job_deltasvm_base_paths = []
job_gkmpredict_base_paths = []
job_snp_sequence_paths = []
for cell_type in cell_types:
    for fold_id in range(1, num_folds+1):
        job_cell_types.append(cell_type)
        job_fold_ids.append(fold_id)
        job_model_base_paths.append(model_base_path)
        job_gkmpredict_base_paths.append(gkmpredict_score_base_path)
        job_deltasvm_base_paths.append(deltasvm_score_base_path)
        job_snp_sequence_paths.append(seq_base_path)
        job_broad_clusters.append(broad_clusters[cell_types.index(cell_type)])

job_params = []
for idx, _ in enumerate(job_cell_types):
    job_params.append({
        "cell_type": job_cell_types[idx],
        "fold_id": job_fold_ids[idx],
        "model_base_path": job_model_base_paths[idx],
        "gkmsvm_base_path": gkmsvm_base_path,
        "gkmpredict_base_path": job_gkmpredict_base_paths[idx],
        "deltasvm_output_base_path": job_deltasvm_base_paths[idx],
        "snp_sequence_base_path": job_snp_sequence_paths[idx],
        "kmer_fasta_path": kmer_fasta_path,
        "broad_cluster": job_broad_clusters[idx]
    })

# # # DEBUG
# job_params = [job_params[0]]