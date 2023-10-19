#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=23:00:00
#SBATCH --mem=500G
#SBATCH -o vireo_demultiplex_snATAC_no_doublet.log

module load python/3.10.2

source ~/python_envs/ENV_genotype_data_packages/bin/activate

vireo -c ../../snATAC_bam_genotyping_cellSNP \
      -d /project/rrg-gturecki/Turecki_Single_Cell_Subject_Genotyping_20220426/newest_version/vcf_hg38_conversion/plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_header_fixed_numeric_chr_names_hg38.vcf.gz \
      -t GT \
      -o /scratch/rdenn/snATAC_vireo_analysis_output_no_doublet \
      --noDoublet \
      --randSeed=1
