#!/bin/bash
#SBATCH --array=1-42
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --time=0-10:00
#SBATCH --account=def-cnagy
#SBATCH -N 1                                           #make output folder for er
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/shells/s_%A_%a.out
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/shells/s_output/sex_%A_%a.er
#SBATCH --mail-user=anjali.chawla@mail.mcgill.ca
#SBATCH --mail-type=ALL

echo "ARR task#: $SLURM_ARRAY_TASK_ID"
#####load variables to use######
DIR_MERGE="/home/anjali5/scratch/graham_merge_analysis/R_scripts"
DIR_ATAC="/home/anjali5/scratch/graham_atac_analysis"
DIR_SHELL="/home/anjali5/scratch/graham_merge_analysis/shells"

cd /home/anjali5/scratch/graham_merge_analysis/shells/
#csv_filename=$( awk "NR==$SLURM_ARRAY_TASK_ID" input_file_list.txt )
line_N=$( awk "NR==$SLURM_ARRAY_TASK_ID"  serial_subjects.txt )

#echo "ARR task#: $SLURM_ARRAY_TASK_ID"
#echo "task: $line_N"

module load StdEnv/2020 r/4.1.0

####Make important Directories
mkdir $DIR_ATAC/$line_N/summary
cp $DIR_ATAC/$line_N/*.tsv.gz $DIR_ATAC/$line_N/summary/
cp $DIR_ATAC/$line_N/*.tsv.gz.tbi $DIR_ATAC/$line_N/summary/
mkdir $DIR_ATAC/$line_N/peaks
mkdir $DIR_ATAC/$line_N/peaks/MACS2

####Call Peaks on fragment files for Signac obj##
Rscript $DIR_MERGE/signac_peaks_pro_metadata_part1.R $line_N

####Generate Raw Matrix & QC barcodes with scATAC-pro##
scATAC-pro -s get_mtx -i $DIR_ATAC/$line_N/summary/fragments.tsv.gz,$DIR_ATAC/$line_N/peaks/MACS2/peaks.bed -c /home/anjali5/scATAC-pro/configure_user.txt -o $DIR_ATAC/$line_N
scATAC-pro -s qc_per_barcode -i $DIR_ATAC/$line_N/summary/fragments.tsv.gz,$DIR_ATAC/$line_N/peaks/MACS2/peaks.bed -c /home/anjali5/scATAC-pro/configure_user.txt -o $DIR_ATAC/$line_N

####Create atac.rds by subsetting cells passed QC in qc_per_barcode (filter less than 5 total unique framgents)##
Rscript $DIR_MERGE/signac_peaks_pro_metadata_part2.R $line_N

####Get final cells that pass QC from ArchR## (has a loop inside for samples #Fix)
Rscript $DIR_MERGE/ArchR_get_barcodes.R

####Prepare for adding amulet multiplets by changing is__cell_barcode in singlecell.csv to 1 for all archR_barcodes##
Rscript $DIR_MERGE/amulet_prep.R $line_N 

####Run Amulet##
bash $DIR_SHELL/amulet.sh

####Add Amulet Multiplet info to atac.rds & subset it fro ArchR##
Rscript $DIR_MERGE/signac_peaks_pro_metadata_part3.R $line_N

####Run CellSNP to demultiplex libraries##
bash $DIR_ATAC/cellsnp.sh 

####Add metadata to ArchR project## (has a function inside for samples #Fix for workflow)
Rscript $DIR_MERGE/Add_metadata_to_ArchR.R 

echo done


