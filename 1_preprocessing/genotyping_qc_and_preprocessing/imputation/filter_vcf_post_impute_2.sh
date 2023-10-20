#!/bin/bash
module load bcftools/1.16
module load picard/2.26.3
module load vcftools/0.1.16

#remove chromosome jumps to sex chromosomes during liftover
grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]' /project/rrg-gturecki/Turecki_Single_Cell_Subject_Genotyping_20220426/newest_version_round_1_only/IMPUTE2_imputation/hg38_post_imputation_liftover_newest/post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_norm.vcf > post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_jump_filter_1.vcf

#change vcf header
bcftools reheader -h new_chr_order_header.txt post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_jump_filter_1.vcf > post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_reorder_header.vcf

#picard reorder vcf based on header
java -jar $EBROOTPICARD/picard.jar SortVcf I=post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_reorder_header.vcf \
                                   O=post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_reordered_fin.vcf

#refilter at INFO 0.8 instead of INFO 0.1
bcftools view -i 'INFO>0.8' post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_reordered_fin.vcf > post_imputation_variants_unified_ref_hg38_info_0.8_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_reordered_fin.vcf

#make a version with MAF 0.5 instead of MAF 0.1
vcftools --vcf post_imputation_variants_unified_ref_hg38_info_0.8_MAF_0.01_hwe_0.000001_sorted_chr_with_info_chr_reordered_fin.vcf --maf 0.05 --recode --recode-INFO-all --out post_imputation_variants_unified_ref_hg38_info_0.8_MAF_0.05_hwe_0.000001_sorted_chr_with_info_chr_reordered_fin.vcf
