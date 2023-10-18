# Calculate LD scores

## Make annotation files
```
for j in $(seq 1 22); do for i in ../Format_Peaks/hg19/*.bed; do python ~/utils/ldsc/make_annot.py --bed-file $i --bimfile ../1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j\.bim --annot-file $(echo $i | cut -f 4 -d '/' | sed "s/\.bed//g").$j\.annot.gz; done; done
```

## Calculate LD scores
```
for j in $(seq 1 22); do for i in ../Format_Peaks/hg19/*.bed; do python ~/utils/ldsc/ldsc.py --print-snps ../1KG/hapmap3_snps/hm.$j\.snp --ld-wind-cm 1.0 --out $(echo $i | cut -f 4 -d '/' | sed "s/\.bed//g").$j --bfile ../1KG/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j --thin-annot --annot $(echo $i | cut -f 4 -d '/' | sed "s/\.bed//g").$j\.annot.gz --l2; done ; done
```
