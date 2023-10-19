#!/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#do linkage disequilibrium pruning

mkdir ld_pruning

plink --ped ../PLINK_starting_files_from_gs_final_report/plus_minus_report.ped \
      --map ../PLINK_starting_files_from_gs_final_report/plus_minus_report.map \
      --geno 0.05 \
      --hwe 0.000001 \
      --indep 50 5 2 \
      --out ld_pruning/ld_pruned

#do relatedness check

plink --ped ../PLINK_starting_files_from_gs_final_report/plus_minus_report.ped \
      --map ../PLINK_starting_files_from_gs_final_report/plus_minus_report.map \
      --genome \
      --extract ld_pruning/ld_pruned.prune.in \
      --out plink_relatedness

sort -nrk 10,10 plink_relatedness.genome > plink.relatedness.genome.sorted
