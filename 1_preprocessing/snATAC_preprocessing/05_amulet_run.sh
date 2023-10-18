#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-cnagy
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/shells/amulet.er
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/shells/amulet.out
#SBATCH --mem=50g
#SBATCH --mail-user=anjali.chawla@mail.mcgill.ca
#SBATCH --mail-type=ALL
DIR="/home/anjali5/scratch/graham_atac_analysis"
module load StdEnv/2020 python/3.8.2 scipy-stack/2020a

#work with fragment file!!!
# A1B1 A2B2 A3B3 A4B4 A5B5 A8B8
for value in A1B1 A2B2 A3B3 A4B4 A5B5 A6B6 A7B7 A8B8 A9B9 A10B10 A11B11 A12B12 A13B13 A14B14 A15B15 A16B16 A17B17 A18B18 A19B19 A20B20 A21B21 A22B22 A23B23 A24B24 A25B25 A26B26 A27B27 A28B28 A29B29 A30B30 A31B31 A32B32 A33B33 A34B34 A35B35 A36B36 A37B37 A38B38 A39B39 A40B40 A41B41 A42B42
do
rm -r ${DIR}/atac_${value}/amulet_archR
mkdir ${DIR}/atac_${value}/amulet_archR
PathToFragments="${DIR}/atac_${value}/fragments.tsv.gz"
PathToBarcode="/home/anjali5/scratch/graham_atac_analysis/atac_${value}/singlecell_archR.csv"
PathToHumanChromosome="/home/anjali5/AMULET/human_autosomes.txt"
PathToRepetetiveElements="/home/anjali5/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed"
OutputDIR="${DIR}/atac_${value}/amulet_archR"
PathToShell="/home/anjali5/AMULET" 

$PathToShell/AMULET.sh $PathToFragments $PathToBarcode $PathToHumanChromosome $PathToRepetetiveElements $OutputDIR $PathToShell 
echo $value
done


