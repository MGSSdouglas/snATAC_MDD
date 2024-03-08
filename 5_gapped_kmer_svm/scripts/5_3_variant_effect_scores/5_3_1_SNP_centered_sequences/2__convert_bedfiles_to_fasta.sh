#!/bin/bash

: <<'COMMENT'
This script converts bedfiles for each SNP to fasta files with respect to hg38.p13 reference genome.
COMMENT

# environment variables set with init_workspace.sh script located in workspace folder
# $GKMSVM_WORKSPACE_DIR
# $GKMSVM_RAW_DATA_DIR
# $GKMSVM_PREPARED_DATA_DIR
# $GKMSVM_MODEL_DIR
# $GKMSVM_TMP_DIR
# $GKMSVM_BIN_DIR

module load bedtools/2.30.0

# get genome fasta file
data_id="hg38/grch38_p13"
data_dir="$GKMSVM_RAW_DATA_DIR/$data_id"
genome_fa="$data_dir/hg38.p13.fa"

# get cdSNPs basepath
cdsnps_basepath="$GKMSVM_PREPARED_DATA_DIR/cdSNPs"

# convert bedfiles for 51bp and 201bp regions 
# centered at cdSNPs per cluster to fasta files
for cluster in ${cdsnps_basepath}/*/ ; do
    for seqlen in 51 201; do
        input_filepath="$cdsnps_basepath/$cluster/${seqlen}bp_seqs/cdsnps.unique.${seqlen}bp.bed"
        output_filepath="$cdsnps_basepath/$cluster/${seqlen}bp_seqs/cdsnps.unique.${seqlen}bp.ref.fa"
        bedtools getfasta -fi "${genome_fa}" -bed "${input_filepath}" -fo "${output_filepath}"
    done
done
    
