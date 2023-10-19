#!/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#check that sex determined from X chr inbreeding coefficient match sex listed in metadata
plink --ped ../PLINK_starting_files_from_gs_final_report/plus_minus_report.ped \
      --map ../PLINK_starting_files_from_gs_final_report/plus_minus_report.map \
      --check-sex \
      --out sex_qc
