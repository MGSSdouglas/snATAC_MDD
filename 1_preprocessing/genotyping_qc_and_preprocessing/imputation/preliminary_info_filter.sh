#!/bin/bash

for chr in $(seq 1 22); do
	mkdir IMPUTE2_results_info_0.1/chr$chr/
	for info_file in IMPUTE2_results/chr$chr/imputed*/*_info; do
		tail -n+2 $info_file | awk '$7 >= 0.1 {print NR}' > gen_lines_tmp.txt
		awk '$7 >= 0.1 {print NR}' $info_file > info_lines_tmp.txt
		awk 'NR == FNR{a[$0]; next};FNR in a' gen_lines_tmp.txt ${info_file/_info/} > IMPUTE2_results_info_0.1/chr$chr/$(basename ${info_file/_info/})
		rm gen_lines_tmp.txt
		awk 'NR == FNR{a[$0]; next};FNR in a' info_lines_tmp.txt $info_file > IMPUTE2_results_info_0.1/chr$chr/$(basename "$info_file")
		rm info_lines_tmp.txt
		head -2 ${info_file/_info/_samples} > IMPUTE2_results_info_0.1/chr$chr/$(basename ${info_file/_info/_samples})
		tail -n+3 ${info_file/_info/_samples} | awk '{$7=$7-1; print $0}' >> IMPUTE2_results_info_0.1/chr$chr/$(basename ${info_file/_info/_samples})
	done
done
