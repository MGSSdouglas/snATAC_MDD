import pdb
import os 
import sys 
import os.path as osp


'''
Configuration
'''

# analysis_type = "subcluster"
analysis_type = "broad"

# cell types to analyze
cell_types = ["Ast1", "Ast2", "Ast3", "Ast4", "End2", "ExN1", \
    "ExN1_L23", "ExN1_L24", "ExN1_L46", "ExN1_L56", "ExN2", "ExN2_L23", \
    "ExN2_L46", "ExN2_L56", "ExN3_L46", "ExN3_L56", "ExN4_L56", "In_LAMP5", \
    "InN3", "In_PV", "In_SST", "In_VIP", "Mic1", "Mic2", "End1", "Oli1", "Oli2", \
    "Oli3", "Oli4", "Oli5", "Oli6", "Oli7", "OPC1", "OPC2", "OPC3", "OPC4"]

broad2sub = {
    "Ast": ["Ast1", "Ast2", "Ast3", "Ast4"],
    "End": ["End2"],
    "ExN": ["ExN1", "ExN1_L23", "ExN1_L24", "ExN1_L46", "ExN1_L56", "ExN2", "ExN2_L23", "ExN2_L46", "ExN2_L56", "ExN3_L46", "ExN3_L56", "ExN4_L56"],
    "InN": ["In_LAMP5", "InN3", "In_PV", "In_SST", "In_VIP"],
    "Mic": ["Mic1", "Mic2", "End1"],
    "Oli": ["Oli1", "Oli2", "Oli3", "Oli4", "Oli5", "Oli6", "Oli7"],
    "OPC": ["OPC1", "OPC2", "OPC3", "OPC4"]
}

broad_clusters = []
for x in cell_types:
    if x in broad2sub["Ast"]:
        broad_clusters.append("Ast")
    elif x in broad2sub["End"]:
        broad_clusters.append("End")
    elif x in broad2sub["ExN"]:
        broad_clusters.append("ExN")
    elif x in broad2sub["InN"]:
        broad_clusters.append("InN")
    elif x in broad2sub["Mic"]:
        broad_clusters.append("Mic")
    elif x in broad2sub["Oli"]:
        broad_clusters.append("Oli")
    elif x in broad2sub["OPC"]:
        broad_clusters.append("OPC")
    else:
        raise ValueError("Cell type {} not found in broad2sub dictionary".format(x))

# cross validation fold count
num_folds = 5

# gapped-kmer svm trained model path
pos_example_type = "fine.before_IPR.60K_thresholding"
model_seqlen = "1001bp"
neg_example_type = "npr1.GC_matched.chrom_matched"
model_train_dataset = "__".join([pos_example_type, model_seqlen, neg_example_type])
model_base_path = f"/home/dcakma3/scratch/mdd-genome2ca/models/{pos_example_type}/{model_seqlen}/{neg_example_type}/cv_gkmsvm"

# path to gapped kmer svm scripts 
gkmsvm_base_path = "/home/dcakma3/projects/def-cnagy/dcakma3/mdd-scripts.v2/mdd-score_snps/lsgkm_kundaje/bin"

'''
I/O file paths
'''

# SNP sequence path
seq_base_path = f"/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/{analysis_type}"

# gkmexplain target path
gkmexplain_score_base_path = f"/home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain"


'''
OLD paths
Example splitted file path:
    /home/dcakma3/scratch/mdd-score_snp_effects/howard_et_al+finemap.TOP/processed_data/union_peaks/51bp/chunks/a1/original/splits
Example gkmexplain score file path of a given split:
    /home/dcakma3/scratch/mdd-score_snp_effects/howard_et_al+finemap.TOP/gkmexplain/a1/original/Ast1/fold_1/imp_score
    /home/dcakma3/scratch/mdd-score_snp_effects/howard_et_al+finemap.TOP/gkmexplain/a1/original/Ast1/fold_1/hyp_imp_score
'''

'''
CURRENT path
Example splitted file path:
    /home/dcakma3/scratch/mdd-prepare_snps_in_peaks/201bp/chunks/a1/original/splits
Example gkmexplain score file path of a given split for Ast1:
    (for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
    (for merged impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original 
    (for splitted hypimpscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/hyp_imp_score/a1/original/splitted_results
    (for merged hypimpscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/hyp_imp_score/a1/original
'''

job_cell_types = []
job_broad_clusters = []
job_fold_ids = []
job_model_base_paths = []
job_gkmexplain_base_paths = []
job_snp_sequence_paths = []
for cell_type in cell_types:
    for fold_id in range(1, num_folds+1):
        job_cell_types.append(cell_type)
        job_fold_ids.append(fold_id)
        job_model_base_paths.append(model_base_path)
        job_gkmexplain_base_paths.append(gkmexplain_score_base_path)
        job_snp_sequence_paths.append(seq_base_path)
        job_broad_clusters.append(broad_clusters[cell_types.index(cell_type)])

job_params = []
for idx, _ in enumerate(job_cell_types):
    job_params.append({
        "cell_type": job_cell_types[idx],
        "fold_id": job_fold_ids[idx],
        "model_base_path": job_model_base_paths[idx],
        "gkmsvm_base_path": gkmsvm_base_path,
        "gkmexplain_base_path": job_gkmexplain_base_paths[idx],
        "snp_sequence_base_path": job_snp_sequence_paths[idx],
        "broad_cluster": job_broad_clusters[idx]
    })


