#!/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#add sort and add cM coordinates to temporary files
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--bfile ../info_0.1/post_imputation_variants_unified_ref_hg19 \
	--make-pgen \
	--sort-vars \
	--out post_imputation_variants_unified_ref_hg19_sorted_tmp

/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--pfile post_imputation_variants_unified_ref_hg19_sorted_tmp \
	--make-bed \
	--out post_imputation_variants_unified_ref_hg19_sorted_tmp

plink \
	--bfile post_imputation_variants_unified_ref_hg19_sorted_tmp \
	--cm-map ../../../SHAPEIT_phasing/1000GP_Phase3/genetic_map_chr@_combined_b37.txt \
	--make-bed \
	--out post_imputation_variants_unified_ref_hg19_sorted_tmp_with_cM

#filter out variants with MAF < 0.01 or hwe < 0.000001
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--bfile post_imputation_variants_unified_ref_hg19_sorted_tmp_with_cM \
	--maf 0.01 \
	--hwe 0.000001 \
	--make-bed \
	--out post_imputation_variants_unified_ref_hg19_fully_filtered

#rm post_imputation_variants_unified_ref_hg19_sorted_tmp*
#rm post_imputation_variants_unified_ref_hg19_sorted_tmp_with_cM*

#make traw file for analysis in R
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--bfile post_imputation_variants_unified_ref_hg19_fully_filtered \
	--export Av \
	--out post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001

#make ped and map files
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
        --bfile post_imputation_variants_unified_ref_hg19_fully_filtered \
	--recode ped \
	--out post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001

#make sorted files

/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
        --bfile post_imputation_variants_unified_ref_hg19_fully_filtered \
	--make-pgen \
        --sort-vars \
	--out post_imputation_variants_unified_info_0.1_MAF_0.01_hwe_0.000001_sorted

/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
        --pfile post_imputation_variants_unified_info_0.1_MAF_0.01_hwe_0.000001_sorted \
        --make-bed \
	--ref-from-fa \
        --fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
        --out post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted
