#!/bin/bash
module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.13
module load picard/2.26.3

mkdir intermediate_vcf_files

#convert plink data to vcf format
/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--bfile ../filtered_plink_hg19/plus_minus_report_filtered_no_dup_pre_vireo \
	--ref-from-fa force \
	--fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
	--normalize \
	--recode vcf-iid \
	--out intermediate_vcf_files/plus_minus_report_filtered_no_dup

#change chromosome name conventions to chrN
bcftools annotate --rename-chrs chr_name_conv_num_to_chrN.txt intermediate_vcf_files/plus_minus_report_filtered_no_dup_no_D_I.vcf > intermediate_vcf_files/plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs.vcf

#lift over to hg38
java -Xmx13G -jar $EBROOTPICARD/picard.jar LiftoverVcf \
        I=intermediate_vcf_files/plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs.vcf \
	O=intermediate_vcf_files/plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_hg38.vcf \
	CHAIN=hg19ToHg38.over.chain \
	REJECT=intermediate_vcf_files/plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_hg38_rejected.vcf \
	R=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
	WARN_ON_MISSING_CONTIG=true

#add missing header fields to hg38 vcf
java -jar $EBROOTPICARD/picard.jar FixVcfHeader \
        I=intermediate_vcf_files/plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_hg38.vcf \
	O=plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_header_fixed_hg38.vcf

#change to numeric chromosome name conventions
bcftools annotate --rename-chrs chr_name_conv_chrN_to_num.txt plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_header_fixed_hg38.vcf > plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_header_fixed_numeric_chr_names_hg38.vcf

#compress and index vcf with numeric chr names so SNP names can be copied to cellSNP output
bgzip -f plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_header_fixed_numeric_chr_names_hg38.vcf
bcftools index -t  plus_minus_report_filtered_no_dup_no_D_I_renamed_chrs_header_fixed_numeric_chr_names_hg38.vcf.gz
