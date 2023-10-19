#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=23:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH -o cellSNP_genotype_snATAC_good_subjects_all_snps_job.log

cellsnp-lite -S snATAC_bam_good_subjects.list -O /scratch/rdenn/cellSNP_genotyping_snATACseq_good_subjects_all_snps -p 32 --cellTAG None --UMItag None --minMAF 0.1 --minCOUNT 20
