#!/bin/bash


# environment variables set with init_workspace.sh script located in workspace folder
# $GKMSVM_WORKSPACE_DIR
# $GKMSVM_ENV_DIR
# $GKMSVM_SCRIPTS_DIR
# $GKMSVM_RAW_DATA_DIR
# $GKMSVM_PREPARED_DATA_DIR
# $GKMSVM_MODEL_DIR
# $GKMSVM_TMP_DIR
# $GKMSVM_BIN_DIR

# load virtualenv
source ${GKMSVM_ENV_DIR}/p38/bin/activate

# create temporary storage directory
mkdir -p "${GKMSVM_TMP_DIR}/5_3_2_score_sequences/2_deltaSVM"

# create deltaSVM result dir 
fpath="${GKMSVM_PREPARED_DATA_DIR}/cdSNP_variant_effect_scores/deltaSVM"
mkdir -p $fpath

# generate all possible non-redundant k-mers with {kmer_length} and 
# save it in {out_file} as FASTA format.
script_path="$GKMSVM_BIN_DIR/lsgkm/scripts/nrkmers.py"
kmer_length="11"
out_file="${fpath}/all_${kmer_length}mers.fa"

# run script
python ${script_path} ${kmer_length} ${out_file}

echo "Script finished."



