#!/bin/bash
#SBATCH --account=def-gturecki
#SBATCH --time=20:00:00
#SBATCH --mem=60G
#SBATCH -o convert_to_hg38_newest.log
module load StdEnv/2020
module load gcc/9.3.0
module load bcftools/1.13
module load picard/2.26.3
module load samtools/1.16

#convert plink data to vcf format

mkdir intermediate_vcf_files

/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--bfile ../hg19_imputed_variants_plink/info_0.1_MAF_0.01_hwe_0.000001/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted \
	--ref-from-fa force \
	--fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
	--recode vcf-iid \
	--out intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted

#change chromosome names
bcftools annotate --rename-chrs chr_name_conv_num_to_chrN.txt intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted.vcf  > intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr.vcf

#add info scores to vcf file
bcftools annotate -c CHROM,FROM,-,REF,ALT,INFO -a combined_info_scores.bed.gz intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr.vcf -h <(echo '##INFO=<ID=INFO,Number=1,Type=String,Description="IMPUTE2 imputation score">') > intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info.vcf

#add info scores for cases where --ref-from-fa switched alleles
bcftools annotate -c CHROM,FROM,-,ALT,REF,INFO -a combined_info_scores.bed.gz intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info.vcf -h <(echo '##INFO=<ID=INFO,Number=1,Type=String,Description="IMPUTE2 imputation score">') > intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_2.vcf

#lift over to hg38
java -Xmx13G -jar $EBROOTPICARD/picard.jar LiftoverVcf \
        I=intermediate_vcf_files/post_imputation_variants_unified_ref_hg19_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_2.vcf \
	O=post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info.vcf \
	CHAIN=hg19ToHg38.over.chain \
	REJECT=post_imputation_variants_unified_ref_REJECT_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info.vcf \
	R=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
	WARN_ON_MISSING_CONTIG=true

#create fasta which chrN style names
sed -r 's/>([0-9]{1,3})/>chr\1/g' /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa > Homo_sapiens.GRCh38.chr_names.fa

samtools faidx Homo_sapiens.GRCh38.chr_names.fa

#normalize variants
bcftools norm \
	--fasta-ref Homo_sapiens.GRCh38.chr_names.fa \
	post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info.vcf > post_imputation_variants_unified_ref_hg38_info_0.1_MAF_0.01_hwe_0.000001_sorted_chr_with_info_norm.vcf
