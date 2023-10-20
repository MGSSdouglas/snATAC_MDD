#!/bin/bash

snp_total=0
for chr in $(seq 1 22); do
	for info_file in IMPUTE2_results/chr$chr/imputed*/*_info; do
		curr_snp_count=$(cut -f7 -d' ' $info_file | awk '$1 >= 0.8' | wc -l)	
		snp_total=$(($snp_total + $curr_snp_count))
	done
done
echo $snp_total
