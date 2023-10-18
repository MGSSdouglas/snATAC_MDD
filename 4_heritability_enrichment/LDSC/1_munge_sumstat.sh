# Format GWAS summary statistics

## Data

We have curated openly available GWAS summary statistics for 
* [Howard Major depressive disorder](https://www.nature.com/articles/s41467-018-03819-3)
* [Als Major depressive disorder](https://pubmed.ncbi.nlm.nih.gov/37464041/)
* [Suicide attempt](https://www.medrxiv.org/content/10.1101/2020.12.01.20241281v1)
* [Insomnia](https://www.nature.com/articles/s41588-022-01124-w)
* [Neuroticism](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006476/sumstats_neuroticism_ctg_format.txt.gz)
* [BMI](https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz)
* [Height](https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz)
* [Attention deficit hyperactivity disorder](https://pgc.unc.edu/for-researchers/download-results/)
* [Schizophrenia](https://pgc.unc.edu/for-researchers/download-results/)
* [Bipolar](https://pgc.unc.edu/for-researchers/download-results/)

## Format GWAS summary statistics

In order to perform LDSC, we use the munge function to format GWAS summary statistics.

```
# MDD
python ~/utils/ldsc/munge_sumstats.py --sumstats PGC_UKB_depression_genome-wide.txt --out MDD --snp MarkerName --N-cas 170756 --N-con 329443 --merge-alleles ../../1KG/w_hm3.snplist --a1 A1 --a2 A2 --frq Freq --p P --signed-sumstats LogOR,0
# MDDALS
python ~/utils/ldsc/munge_sumstats.py --sumstats
daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_MVPaf_HRC_MAF01.gz --out MDD_ALS --snp SNP --N-cas 294322 --N-con 741438 --merge-alleles ~/projects/def-cnagy/For_Wenmin/scMDD_GWAS_WZ/1KG/w_hm3.snplist --a1 A1 --a2 A2 --frq FRQ_A_294322 --p P --signed-sumstats OR,1 
# SA
python ~/utils/ldsc/munge_sumstats.py --sumstats daner_model2_062620_eur.neff.qc2.80.gz --daner-n --out SA --a1 A1.y --a2 A2.y
# INS
python ~/utils/ldsc/munge_sumstats.py --sumstats insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz --out INS --N-col NMISS
# NEU
python ~/utils/ldsc/munge_sumstats.py --sumstats sumstats_neuroticism_ctg_format.txt.gz --N-col N --out NEU --merge-alleles ../../1KG/w_hm3.snplist --a1 A1 --a2 A2 --p P --signed-sumstats Z,0 --snp RSID
# BMI
python ~/utils/ldsc/munge_sumstats.py --sumstats Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz --out BMI --a1 Tested_Allele --a2 Other_Allele
# HGT
python ~/utils/ldsc/munge_sumstats.py --sumstats Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz --out HGT --a1 Tested_Allele --a2 Other_Allele
# ADHD 
python ~/utils/ldsc/munge_sumstats.py --sumstats daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz --daner-n --out ADHD
# SCZ
python ~/utils/ldsc/munge_sumstats.py --sumstats PGC3_SCZ_wave3_public.v2.tsv.gz --N-cas 67390 --N-con 94015 --out SCZ --merge-alleles ../../1KG/w_hm3.snplist --a1 A1 --a2 A2 --frq FRQ_A_67390 --p P --signed-sumstats OR,1
# BIP
python ~/utils/ldsc/munge_sumstats.py --sumstats daner_PGC_BIP32b_mds7a_0416a.gz --daner-n --out BIP
```


