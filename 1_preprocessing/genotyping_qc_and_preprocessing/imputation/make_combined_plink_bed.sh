#/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

for chr in $(seq 1 22); do
        for sample_file in ../../IMPUTE2_results_info_0.8/chr$chr/imputed*_samples; do
		/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
			--gen ${sample_file/_samples/} ref-first \
			--sample $sample_file  \
			--oxford-single-chr $chr \
			--make-bed \
			--out chr"$chr"_"$(basename "${sample_file/_samples/}")"
	done
done

ls -1q chr*.bed > bed_file.list
ls -1q chr*.bim > bim_file.list
ls -1q chr*.fam > fam_file.list
paste bed_file.list bim_file.list fam_file.list > bfiles.list
head -n -1 bfiles.list > files_to_merge.list

plink \
	--bfile chr9_imputed_99000000_102000000 \
	--merge-list files_to_merge.list \
	--out post_imputation_variants_unified_hg19

/project/rrg-gturecki/Software_Installations_and_Databases/plink2 \
	--bfile post_imputation_variants_unified_hg19 \
	--ref-from-fa force \
	--fa /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
	--normalize \
	--make-bed \
	--out post_imputation_variants_unified_ref_hg19 
