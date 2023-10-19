#!/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#convert from lgen to ped/map format, while adding cM positions of SNPs and removing indels
plink --lfile plink_lgen_no_cM/plus_minus_report \
      --cm-map ../SHAPEIT_phasing/1000GP_Phase3/genetic_map_chr@_combined_b37.txt \
      --recode \
      --out plus_minus_report
