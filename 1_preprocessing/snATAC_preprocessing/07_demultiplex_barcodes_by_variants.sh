#!/bin/bash
#SBATCH --time=0-20:00:00
#SBATCH --account=def-cnagy
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/shells/cellsnp.er
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/shells/cellsnp.out
#SBATCH --mem=50g
#SBATCH --mail-user=anjali.chawla@mail.mcgill.ca
#SBATCH --mail-type=ALL

#module load nixpkgs/16.09 python/3.7.4
#module load scipy-stack
#module load samtools


conda init bash
source ~/.bashrc
conda activate /home/anjali5/conda3/envs/CSP
DIR_ATAC="/home/anjali5/scratch/graham_atac_analysis"
DIR_MERGE="/home/anjali5/scratch/graham_merge_analysis"
DIR_BAM="/home/anjali5/projects/def-cnagy/Nanuq_Aug21_snATAC/bams"

#########based on refSNPs all chromosomes pileup#########
##0.1,20## noerror is the best
##0.05,10## 1 error
##0.02,5## many error
##0.01,2## many error
##0.05,20##1 error less cells

#module load samtools

for value in A1B1 A2B2 A3B3 A4B4 A5B5 A6B6 A7B7 A8B8 A9B9 A10B10 A11B11 A12B12 A13B13 A14B14 A15B15 A16B16 A17B17 A18B18 A19B19 A20B20 A21B21 A22B22 A23B23 A24B24 A25B25 A26B26 A27B27 A28B28 A29B29 A30B30 A31B31 A32B32 A33B33 A34B34 A35B35 A36B36 A37B37 A38B38 A39B39 A40B40 A41B41 A42B42
do
	awk 'NR>1' $DIR_ATAC/atac_$value/archR_barcodes > $DIR_ATAC/atac_$value/archR_barcodes_for_cellsnp
	cellsnp-lite -s $DIR_BAM/*.$value.bam -b $DIR_ATAC/atac_$value/archR_barcodes_for_cellsnp -R $DIR_MERGE/human_snp/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz -O $DIR_ATAC/atac_$value/snp_to_check --minMAF 0.1 --minCOUNT 20 --UMItag None -p 32 --genotype --gzip --minMAPQ 30 --inclFLAG 2
	vireo -c $DIR_ATAC/atac_$value/snp_to_check -N 2 -o $DIR_ATAC/atac_$value/vireo_to_check
	cd $DIR_ATAC/atac_$value/vireo/
	cut -f 1,2 donor_ids.tsv > dtemp
	grep "donor0" dtemp | cut -f 1 > donor0
	grep "donor1" dtemp | cut -f 1 > donor1
	grep "doublet" dtemp | cut -f 1 > doublet
	cat donor1 donor1 > donors.txt
	echo $value 
	#head donors.txt
	wc -l donors.txt
done
conda deactivate 


#GTbarcode -i /home/anjali5/scratch/graham_atac_analysis/snp_to_check/cellSNP.cells.vcf.gz -o /home/anjali5/scratch/graham_atac_analysis/atac_bam/GT_barcodes.tsv --randSeed 1



