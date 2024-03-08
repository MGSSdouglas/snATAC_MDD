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
gkmexplain_score_base_path=$4
seq_base_path=$5
gkmexplain="$6/gkmexplain"
clustering_type=$7


####################################
# 201bp, a1, original
####################################
seq_len="201bp"
seq_type="original"
snp_type="a1"

split_file_base_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/chunks/${snp_type}/${seq_type}/splits"
#/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/broad/Ast/Ast1/201bp/chunks/a1/original/splits

# generate importance score using 16 cores
score_type="imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-00.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-01.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-02.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-03.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-04.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-05.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-06.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-07.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-08.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-09.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-10.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-11.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-12.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-13.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-14.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain

echo "Done: $seq_len $snp_type $seq_type $score_type"

# generate hypothetical importance score using 16 cores
score_type="hyp_imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-00.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-01.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-02.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-03.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-04.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-05.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-06.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-07.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-08.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-09.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-10.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-11.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-12.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-13.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-14.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain

echo "Done with $snp_type $seq_type $score_type"

####################################
# 201bp, a1, shuffled
####################################
seq_len="201bp"
seq_type="shuffled"
snp_type="a1"

split_file_base_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/chunks/${snp_type}/${seq_type}/splits"
#/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/broad/Ast/Ast1/201bp/chunks/a1/original/splits

# generate importance score using 16 cores
score_type="imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain

echo "Done: $seq_len $snp_type $seq_type $score_type"

# generate hypothetical importance score using 16 cores
score_type="hyp_imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain

echo "Done with $snp_type $seq_type $score_type"

####################################
# 201bp, a2, original
####################################
seq_len="201bp"
seq_type="original"
snp_type="a2"

split_file_base_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/chunks/${snp_type}/${seq_type}/splits"
#/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/broad/Ast/Ast1/201bp/chunks/a1/original/splits

# generate importance score using 16 cores
score_type="imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-00.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-01.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-02.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-03.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-04.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-05.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-06.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-07.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-08.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-09.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-10.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-11.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-12.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-13.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-14.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain

echo "Done: $seq_len $snp_type $seq_type $score_type"

# generate hypothetical importance score using 16 cores
score_type="hyp_imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-00.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-01.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-02.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-03.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-04.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-05.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-06.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-07.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-08.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-09.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-10.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-11.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-12.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-13.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-14.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}-*.gkmexplain

echo "Done with $snp_type $seq_type $score_type"

####################################
# 201bp, a2, shuffled
####################################
seq_len="201bp"
seq_type="shuffled"
snp_type="a2"

split_file_base_path="${seq_base_path}/${clustering_type}/${cell_type}/${seq_len}/chunks/${snp_type}/${seq_type}/splits"
#/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/broad/Ast/Ast1/201bp/chunks/a1/original/splits

# generate importance score using 16 cores
score_type="imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.gkmexplain" & \
${gkmexplain} "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain

echo "Done: $seq_len $snp_type $seq_type $score_type"

# generate hypothetical importance score using 16 cores
score_type="hyp_imp_score"
#(for splitted impscore results) /home/dcakma3/scratch/mdd-score_snps.v2/{analysis_type}/gkmexplain/Ast/Ast1/fold_1/imp_score/a1/original/splitted_results 
outfile_base_path="${gkmexplain_score_base_path}/${clustering_type}/${cell_type}/fold_${fold_id}/${seq_len}/${snp_type}/${seq_type}"
mkdir -p "${outfile_base_path}/${score_type}"
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-00.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-01.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-02.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-03.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-04.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-05.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-06.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-07.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-08.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-09.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-10.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-11.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-12.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-13.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-14.gkmexplain" & \
${gkmexplain} -m 1  "$split_file_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.fa" ${trained_model_path} \
    "$outfile_base_path/$score_type/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-15.gkmexplain" & \
wait
merged_file_path="$outfile_base_path/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}.${score_type}.gkmexplain"
for chunk in ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain
do
    cat $chunk >> $merged_file_path
done
rm ${outfile_base_path}/${score_type}/snps.hg38_information.${seq_len}.${snp_type}.${seq_type}-*.gkmexplain

echo "Done with $snp_type $seq_type $score_type"