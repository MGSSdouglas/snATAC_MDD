#!/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#filter plink data for only subjects passing QC and with MAF > 0.01, call rate > 95%, hwe p-value > 0.000001
mkdir plink_filtered_with_duplicates_pre_vireo

/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
      --ped ../PLINK_starting_files_from_gs_final_report/plus_minus_report.ped \
      --map ../PLINK_starting_files_from_gs_final_report/plus_minus_report.map \
      --remove bad_subjects_for_remove_pre_vireo_all.txt \
      --chr 1-22 \
      --maf 0.01 \
      --geno 0.05 \
      --hwe 0.000001 \
      --ref-from-fa force \
      --fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
      --normalize \
      --make-bed \
      --out plink_filtered_with_duplicates_pre_vireo/plus_minus_report_filtered_with_dup_pre_vireo

#get list of duplicated SNPs based on position
plink --bfile plink_filtered_with_duplicates_pre_vireo/plus_minus_report_filtered_with_dup_pre_vireo \
      --list-duplicate-vars ids-only suppress-first \
      --out plink_duplicates_pre_vireo_filters

#remove duplicates from dataset
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
      --bfile plink_filtered_with_duplicates_pre_vireo/plus_minus_report_filtered_with_dup_pre_vireo \
      --exclude plink_duplicate_pre_vireo_filters.dupvar \
      --ref-from-fa force \
      --fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
      --normalize \
      --recode ped \
      --out plus_minus_report_filtered_no_dup_pre_vireo

#make bed format output as well
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
      --ped plus_minus_report_filtered_no_dup_pre_vireo.ped \
      --map plus_minus_report_filtered_no_dup_pre_vireo.map \
      --make-bed \
      --out plus_minus_report_filtered_no_dup_pre_vireo
