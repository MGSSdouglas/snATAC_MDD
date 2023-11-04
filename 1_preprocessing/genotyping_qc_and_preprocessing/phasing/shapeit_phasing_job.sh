#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=50:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH -o shapeit_phasing_job.log
module load shapeit/2.r904

#make directories for output and logs
mkdir shapeit_phased_data
mkdir shapeit_logs

#iterate over chromosomes
for chr in $(seq 1 22); do

	#define input files
	curr_bed=filtered_sep_by_chr_hg19/chr$chr/filtered_plink_chr$chr.bed
	curr_bim=filtered_sep_by_chr_hg19/chr$chr/filtered_plink_chr$chr.bim
	curr_fam=filtered_sep_by_chr_hg19/chr$chr/filtered_plink_chr$chr.fam

	#make output directory
	mkdir shapeit_phased_data/chr$chr
	
	#check for SNPs in study but not in 1000G or strand issues
	shapeit -check \
		--input-bed $curr_bed $curr_bim $curr_fam \
                --input-map 1000GP_Phase3/genetic_map_chr"$chr"_combined_b37.txt \
                --input-ref 1000GP_Phase3/1000GP_Phase3_chr$chr.hap.gz 1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz  1000GP_Phase3/1000GP_Phase3.sample \
		--output-log shapeit_logs/shapeit_check_chr$chr

	#run shapeit with problematic SNPs removed
	shapeit --input-bed $curr_bed $curr_bim $curr_fam \
		--input-map 1000GP_Phase3/genetic_map_chr"$chr"_combined_b37.txt \
		--input-ref 1000GP_Phase3/1000GP_Phase3_chr$chr.hap.gz 1000GP_Phase3/1000GP_Phase3_chr$chr.legend.gz  1000GP_Phase3/1000GP_Phase3.sample \
		--output-max shapeit_phased_data/chr$chr/phased_chr$chr.haps shapeit_phased_data/chr$chr/phased_chr$chr.sample \
		--exclude-snp shapeit_logs/shapeit_check_chr$chr.snp.strand.exclude \
		--states 200 \
		--burn 15 \
		--main 40 \
		--output-log shapeit_logs/shapeit_chr$chr.log \
		--seed 1 
done
