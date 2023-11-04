#!/bin/bash

mkdir filtered_sep_by_chr_hg19

#separate filtered plink genotyping data by chromosome
for chr in $(seq 1 22); do
	mkdir filtered_sep_by_chr_hg19/chr$chr
	/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
		--bfile ../filtered_plink_hg19/plus_minus_report_filtered_no_dup_post_vireo \
		--ref-from-fa force \
		--fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
		--normalize \
		--chr $chr \
		--make-bed \
		--out filtered_sep_by_chr_hg19/chr$chr/filtered_plink_chr$chr ;
done
