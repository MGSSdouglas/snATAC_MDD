#!/bin/bash
#SBATCH --mem=200G
#SBATCH --account=def-cnagy
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
#SBATCH -N 1
#SBATCH --output=/home/anjali5/scratch/graham_merge_analysis/R_scripts/ldsc0.out
#SBATCH --error=/home/anjali5/scratch/graham_merge_analysis/R_scripts/ldsc0.er
#SBATCH --mail-user=anjali.chawla@mail.mcgill.ca
#SBATCH --mail-type=ALL


source ~/opt/conda/bin/activate
cd ~/projects/def-cnagy/anjali5/ldsc_ldsc/ldsc
conda activate ldsc

cd /home/anjali5/scratch/LDSC_for_DAR
#cell-type list
sed "s/\"//g" combinedDAR.csv | cut -f 4 -d ',' | sed 1d | sort | uniq > celltypes.lst

#peaks per cell-type
mkdir hg38
while read a; do for i in combinedDAR.csv; do sed "s/\"//g" $i | cut -f 3,4 -d ',' | awk -F ',' -v var=$a '$2==var{print}' | cut -f 1 -d ',' | awk -F '-' '{print $1"\t"$2"\t"$3"\t"$0}' > ./hg38/$i\_$a\.bed ; done; done < celltypes.lst

#peaks acorss cell-types
for i in combinedDAR.csv; do cut -f 3 $i -d ',' | sed 1d | sed "s/\"//g" | sort | uniq | sed "s/-/\t/g" > ./hg38/$i\.bed; done

#convert to hg19
mkdir hg19
module load  StdEnv/2020 kentutils/401
for i in ./hg38/*.bed; do liftOver $i ~/scratch/LDSC_for_DAR/hg38ToHg19.over.chain.gz $(echo $i | sed "s/38/19/g") unmapped ; done

#LDScore
mkdir ldscore
cd ~/scratch/LDSC_for_DAR/ldscore
for j in $(seq 1 22); do for i in ~/scratch/LDSC_for_DAR/hg19/*.bed; do python ~/projects/def-cnagy/anjali5/ldsc_ldsc/ldsc/make_annot.py --bed-file $i --bimfile ~/projects/def-cnagy/For_Wenmin/scMDD_GWAS_WZ/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j\.bim --annot-file $(echo $i | cut -f 7 -d '/' | sed "s/\.bed//g").$j\.annot.gz; done; done
for j in $(seq 20 22); do for i in ~/scratch/LDSC_for_DAR/hg19/*.bed; do python ~/projects/def-cnagy/anjali5/ldsc_ldsc/ldsc/ldsc.py --print-snps ~/projects/def-cnagy/For_Wenmin/scMDD_GWAS_WZ/1KG/hapmap3_snps/hm.$j\.snp --ld-wind-cm 1.0 --out $(echo $i | cut -f 7 -d '/' | sed "s/\.bed//g").$j --bfile ~/projects/def-cnagy/For_Wenmin/scMDD_GWAS_WZ/1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j --thin-annot --annot $(echo $i | cut -f 7 -d '/' | sed "s/\.bed//g").$j\.annot.gz --l2; done ; done


#trait enrichment
cd ~/scratch/LDSC_for_DAR
mkdir cts
cd cts
awk '{print$1"\t../ldscore/CombinedDAR.csv_"$1".,../ldscore/CombinedDAR.csv."}' ~/scratch/LDSC_for_DAR/celltypes.lst > combinedDAR.csv.ldcts
for i in ~/scratch/scMDD_GWAS_WZ/Format_GWAS/*.sumstats.gz; do python ~/projects/def-cnagy/anjali5/ldsc_ldsc/ldsc/ldsc.py --h2-cts $i --ref-ld-chr ~/projects/def-cnagy/For_Wenmin/scMDD_GWAS_WZ/1KG/1000G_EUR_Phase3_baseline/baseline. --out $(echo $i | cut -d '/' -f 7 | sed "s/.sumstats.gz//g")_combinedDARs.csv  --ref-ld-chr-cts combinedDARs.csv.ldcts --w-ld-chr ~/projects/def-cnagy/For_Wenmin/scMDD_GWAS_WZ/1KG/weights_hm3_no_hla/weights.; done

