#!/bin/bash
#SBATCH --mem=200G
#SBATCH --account=def-cnagy
#SBATCH --cpus-per-task=16
#SBATCH --time=0-8:00:00
#SBATCH -N 1
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/R_scripts/step2.out
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/R_scripts/step2.er
#SBATCH --mail-user=anjali.chawla@mail.mcgill.cae
#SBATCH --mail-type=ALL

Rscript ~/scratch/graham_merge_analysis/R_scripts/create_snATAC_project.R
