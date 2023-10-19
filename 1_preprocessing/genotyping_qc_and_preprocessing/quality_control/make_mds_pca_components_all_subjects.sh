#!/bin/bash
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

#do ld pruning before pca

mkdir pca_ld_prune_lists_all_subjects

plink --bfile ../filtered_plink_hg19/plus_minus_report_filtered_no_dup_post_vireo_pos_id_no_MHC \
      --indep 50 5 2 \
      --out pca_ld_prune_lists_all_subjects/plus_minus_report_filt

mkdir plink_pca_pruned_filtered_all_subjects

plink --bfile ../filtered_plink_hg19/plus_minus_report_filtered_no_dup_post_vireo_pos_id_no_MHC \
      --extract pca_ld_prune_lists_all_subjects/plus_minus_report_filt.prune.in \
      --make-bed \
      --out plink_pca_pruned_filtered_all_subjects/plus_minus_report_prune_filt_for_pca

#do pca and make mds plot for visualization

mkdir plink_eigenstrat_pca_all_subjects

plink --bfile plink_pca_pruned_filtered_all_subjects/plus_minus_report_prune_filt_for_pca \
      --pca \
      --out plink_eigenstrat_pca_all_subjects/plus_minus_report_prune_filt_for_pca

mkdir plink_mds_all_subjects

plink --bfile plink_pca_pruned_filtered_all_subjects/plus_minus_report_prune_filt_for_pca \
      --cluster \
      --mds-plot 10 \
      --out plink_mds_all_subjects/plus_minus_report_prune_filt_for_pca
