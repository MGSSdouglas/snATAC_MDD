## Fine-mapping MDD GWAS

We curated the MDD GWAS summary statistics based on UK Biobank from Howard et al. 

```
wget https://datashare.ed.ac.uk/download/DS_10283_3203.zip
```

To extract z-scores and match ref/alt alleles with LD information, we formatted GWAS summary statistics:

```
python matchss.py --rss PGC_UKB_depression_genome-wide.txt --prefix MDD --save MDD --idir ukb/idx/ --rsid MarkerName --A1 A1 --A2 A2 --BETA LogOR --SE StdErrLogOR --delim ' '
```

To finemap with SparsePro:

```
cd MDD
for i in $(seq 1 22); do python ../sparsepro_ukb.py --ukb ~/utils/ukb/lst/$i\.lst --zdir MDD_$i\.z --N $(cat MDD.N) --save sparsepro --prefix MDD_$i --verbose --LDdir ~/scratch/UKBBLD/ --K 5 --cthres 0.9 --ethres 50; done
```

The obtained fine-mapping results were saved into sparsepro folder. Two files with suffix `.pip` and `.cs` were generated for each chromosome, where `.pip` file contains variant level information and `.cs` file contains effect level information. 

To annotate each effect group to the closest gene:

```
cd MDD
python ../map2gene.py --gtf ~/utils/gencode.v19.protein_coding.gtf --rs ~/utils/ukb/idx/allrsid.dict --cs sparsepro/MDD_ --save MDD_cgene.txt --eQTL ~/scratch/SharePro_gxe/dat/anno/GTEx.sQTL.rsid.dict --sQTL ~/scratch/SharePro_gxe/dat/anno/GTEx.eQTL.rsid.dict
```

The result was saved to MDD_cgene.txt.
