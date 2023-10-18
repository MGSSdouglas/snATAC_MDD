# Estimate heritability

## Calculate broad cell groups heritability
```
while read a; do for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2 $i --ref-ld-chr ../Old_LD_Score/$a\. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_$a\_broad --w-ld-chr ../1KG/weights_hm3_no_hla/weights. ; done; done < ../Format_Peaks/broad.lst
```

## Calculate subtype heritability
```
while read a; do for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2 $i --ref-ld-chr ../Old_LD_Score/$a\. --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_$a\_sub --w-ld-chr ../1KG/weights_hm3_no_hla/weights. ; done; done < ../Format_Peaks/sub.lst
```

## Calculate trait heritability
```
for i in ../Format_GWAS/*.sumstats.gz; do python ~/utils/ldsc/ldsc.py --h2 $i --ref-ld-chr ../1KG/eur_w_ld_chr/ --out $(echo $i | cut -d '/' -f 3 | sed "s/.sumstats.gz//g")_base --w-ld-chr ../1KG/weights_hm3_no_hla/weights. ; done
```
