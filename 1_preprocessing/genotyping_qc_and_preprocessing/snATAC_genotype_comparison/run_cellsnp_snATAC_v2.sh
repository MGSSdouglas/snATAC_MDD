#!/bin/bash
#SBATCH --account=def-cnagy
#SBATCH --time=23:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH -o cellSNP_genotype_snATAC_good_subjects_all_snps_job_v2.log

export PATH="/project/rrg-gturecki/Software_Installations_and_Databases/cellsnp-lite:$PATH"

cellsnp-lite -S snATAC_bam_good_subjects.list -O cellSNP_genotyping_snATACseq_good_subjects_all_snps -p 32 --cellTAG None --UMItag None --minMAF 0.1 --minCOUNT 20 --genotype
