#!/bin/bash
#SBATCH --account=ctb-liyue
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=01:20:00
#SBATCH --output=./log/%x.out
#SBATCH --mail-user=doruk.cakmakci@mail.mcgill.ca
#SBATCH --mail-type=FAIL

echo "Job started at: `date`"
echo "Job ID: $SLURM_JOBID"

# input parameters
cell_type=$1
clustering_type=$9
fold_id=$2
trained_model_path="$3/${cell_type}/fold_${fold_id}/trained.model.txt"
deltasvm_score_base_path=$4
gkmpredict_score_base_path=$5
seq_base_path=$6
gkmexplain="$7/gkmexplain"
gkmpredict="$7/gkmpredict"

# define input path
kmer_fasta_path=$8

# define output path
kmer_score_path="${deltasvm_score_base_path}/${clustering_type}/${cell_type}/kmer_scores/${cell_type}/fold_${fold_id}"
mkdir -p $kmer_score_path
kmer_score_path="${kmer_score_path}/kmer_scores.txt"

# score nkmers using the model
$gkmpredict -T 16 $kmer_fasta_path $trained_model_path $kmer_score_path
