#!/bin/bash
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16
#SBATCH --time=0-10:00
#SBATCH --account=def-cnagy
#SBATCH -N 1
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/R_scripts/rscript.out
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/R_scripts/rscript.er
#SBATCH --mail-user=anjali.chawla@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09 gcc/7.3.0
module load r/3.6.1
module load StdEnv/2020 r/4.1.0
module load mugqic/homer/4.11

PATH=$PATH:/home/anjali5/.local/bin/bin:$PATH #install homer locally

folder="/home/anjali5/projects/def-gturecki/anjali5/Finalized_DAR/SubClusters_Motifs_Homer"
bed=".bed"
ext=".out"
input="/home/anjali5/projects/def-gturecki/anjali5/Finalized_DAR/SubClusters_Motifs_Homer"
cd $input
for filename in *.bed;
do
echo $filename
output="$filename$ext"
#cur="$filename$bed"
findMotifsGenome.pl $filename hg38 $folder/$output -p 32 
#annotatePeaks.pl $filename hg38 > $annot/$filename.txt   #annotate peaks
echo $filename$bed
echo $output
done
