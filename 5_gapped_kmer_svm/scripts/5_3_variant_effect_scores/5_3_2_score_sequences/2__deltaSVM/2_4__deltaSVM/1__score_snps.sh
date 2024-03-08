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
fold_id=$2
trained_model_path="$3/${cell_type}/fold_${fold_id}/trained.model.txt"
deltasvm_score_base_path=$4
gkmpredict_score_base_path=$5
seq_base_path=$6
gkmexplain="$7/gkmexplain"
gkmpredict="$7/gkmpredict"
deltasvm_script_path=$8
clustering_type=$9

outfile_base="${deltasvm_score_base_path}/${clustering_type}/${cell_type}/deltaSVM_scores"


####################################
# 201bp, original
####################################
seq_len="201bp"
seq_type="original"

# define input sequence path for effect (A1) allele
snp_type="a1"
effect_allele_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/snps.hg38_information.${seq_len}.${snp_type}"
if [ "$seq_type" == "original" ]; then
    effect_allele_path="${effect_allele_path}.no_lb.fa"
elif [ "$seq_type" == "shuffled" ]; then
    effect_allele_path="${effect_allele_path}.shuffled.no_lb.fa"
fi

# define input sequence path for non-effect (A2) allele
snp_type="a2"
non_effect_allele_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/snps.hg38_information.${seq_len}.${snp_type}"
if [ "$seq_type" == "original" ]; then
    non_effect_allele_path="${non_effect_allele_path}.no_lb.fa"
elif [ "$seq_type" == "shuffled" ]; then
    non_effect_allele_path="${non_effect_allele_path}.shuffled.no_lb.fa"
fi

# define kmer scores path 
kmer_scores_path="${deltasvm_score_base_path}/${clustering_type}/${cell_type}/kmer_scores/${cell_type}/fold_${fold_id}/kmer_scores.txt"

# define output file path for deltaSVM scores
out_file_path="${outfile_base}/${cell_type}/fold_${fold_id}/"
mkdir -p $out_file_path
out_file_path="${out_file_path}/snps.hg38_information.${seq_len}.${seq_type}.a1_a2.deltasvm"

# run deltaSVM
perl $deltasvm_script_path $non_effect_allele_path $effect_allele_path $kmer_scores_path $out_file_path


####################################
# 201bp, shuffled
####################################
seq_len="201bp"
seq_type="shuffled"

# define input sequence path for effect (A1) allele
snp_type="a1"
effect_allele_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/snps.hg38_information.${seq_len}.${snp_type}"
if [ "$seq_type" == "original" ]; then
    effect_allele_path="${effect_allele_path}.no_lb.fa"
elif [ "$seq_type" == "shuffled" ]; then
    effect_allele_path="${effect_allele_path}.shuffled.no_lb.fa"
fi

# define input sequence path for non-effect (A2) allele
snp_type="a2"
non_effect_allele_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/snps.hg38_information.${seq_len}.${snp_type}"
if [ "$seq_type" == "original" ]; then
    non_effect_allele_path="${non_effect_allele_path}.no_lb.fa"
elif [ "$seq_type" == "shuffled" ]; then
    non_effect_allele_path="${non_effect_allele_path}.shuffled.no_lb.fa"
fi

# define kmer scores path 
kmer_scores_path="${deltasvm_score_base_path}/${clustering_type}/${cell_type}/kmer_scores/${cell_type}/fold_${fold_id}/kmer_scores.txt"

# define output file path for deltaSVM scores
out_file_path="${outfile_base}/${cell_type}/fold_${fold_id}/"
mkdir -p $out_file_path
out_file_path="${out_file_path}/snps.hg38_information.${seq_len}.${seq_type}.a1_a2.deltasvm"

# run deltaSVM
perl $deltasvm_script_path $non_effect_allele_path $effect_allele_path $kmer_scores_path $out_file_path