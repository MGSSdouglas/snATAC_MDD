#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=01:20:00
#SBATCH --output=./log/%x.out
#SBATCH --mail-type=FAIL

echo "Job started at: `date`"
echo "Job ID: $SLURM_JOBID"

# input parameters
cluster=$1
fold_id=$2
trained_model_path=$3
gkmpredict="$4/gkmpredict"
cdsnp_dir=$5
cdsnp_score_dir=$6


####################################
# 201bp, a1, original
####################################
seq_len="201bp"
seq_type="original"
snp_type="a1"

# define input path 
input_path="${cdsnp_basepath}/${cluster}/${seq_len}_seqs/cdsnps.unique.${seq_len}.${snp_type}.${seq_type}.fa"

# define output filepath
outfile_base_path="$cdsnp_score_dir/${cluster}/${fold_id}/${seq_len}/${seq_type}/"
mkdir -p $outfile_base_path
outfile_path="$outfile_base_path/cdsnps.unique.${seq_len}.${seq_type}.${snp_type}.gkmpredict"

# run gkmpredict
$gkmpredict -T 16 $input_path $trained_model_path $outfile_path




####################################
# 201bp, a1, shuffled
####################################
seq_len="201bp"
seq_type="shuffled"
snp_type="a1"

# define input path 
input_path="${cdsnp_basepath}/${cluster}/${seq_len}_seqs/cdsnps.unique.${seq_len}.${snp_type}.${seq_type}.fa"

# define output filepath
outfile_base_path="$cdsnp_score_dir/${cluster}/${fold_id}/${seq_len}/${seq_type}/"
mkdir -p $outfile_base_path
outfile_path="$outfile_base_path/cdsnps.unique.${seq_len}.${seq_type}.${snp_type}.gkmpredict"

# run gkmpredict
$gkmpredict -T 16 $input_path $trained_model_path $outfile_path




####################################
# 201bp, a2, original
####################################
seq_len="201bp"
seq_type="original"
snp_type="a2"

# define input path 
input_path="${cdsnp_basepath}/${cluster}/${seq_len}_seqs/cdsnps.unique.${seq_len}.${snp_type}.${seq_type}.fa"

# define output filepath
outfile_base_path="$cdsnp_score_dir/${cluster}/${fold_id}/${seq_len}/${seq_type}/"
mkdir -p $outfile_base_path
outfile_path="$outfile_base_path/cdsnps.unique.${seq_len}.${seq_type}.${snp_type}.gkmpredict"

# run gkmpredict
$gkmpredict -T 16 $input_path $trained_model_path $outfile_path




####################################
# 201bp, a2, shuffled
####################################
seq_len="201bp"
seq_type="shuffled"
snp_type="a2"

# define input path 
input_path="${cdsnp_basepath}/${cluster}/${seq_len}_seqs/cdsnps.unique.${seq_len}.${snp_type}.${seq_type}.fa"

# define output filepath
outfile_base_path="$cdsnp_score_dir/${cluster}/${fold_id}/${seq_len}/${seq_type}/"
mkdir -p $outfile_base_path
outfile_path="$outfile_base_path/cdsnps.unique.${seq_len}.${seq_type}.${snp_type}.gkmpredict"

# run gkmpredict
$gkmpredict -T 16 $input_path $trained_model_path $outfile_path