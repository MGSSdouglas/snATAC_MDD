#!/bin/bash
#SBATCH --mem=200G
#SBATCH --account=def-cnagy
#SBATCH --cpus-per-task=16
#SBATCH --time=0-10:00:00
#SBATCH -N 1
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/R_scripts/step3.out
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/R_scripts/step3.er
#SBATCH --mail-user=anjali.chawla@mail.mcgill.cae
#SBATCH --mail-type=ALL

Rscript ~/scratch/graham_merge_analysis/R_scripts/compute_silhouette_for_clustering.R
