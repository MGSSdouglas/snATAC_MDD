## Main 
This directory contains scripts for the multi-modal visualization of MDD SNPs. Please note we adapted ArchR scripts (Granja et. al, 2021) and added several functionalities to 1. add external data (such as, histone peaks, 3D chromatin loops, external peak-to-gene linkages). 2. highlight SNPs and target genes. 3. add SNPs in LD with signficant MDD SNPs. 4. easy plotting and coloring of cell-type or cluster-specific peaks.

## To prepare DLPFC histone modification peaks (https://www.synapse.org/#!Synapse:syn12245061/tables)

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/6_snp_visualization/4_per_subject_histone_for_snp_tracks.R)

## To prepare DLPFC 3D HiC loops (DLPFC Hi-C 3D chromatin data; PMID: 34172755)

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/6_snp_visualization/3_getHiC.R)

## To prepare peak-to-peak and peak-to-gene inputs

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/6_snp_visualization/2_get_input_p2p_p2g.R)

## To prepare for multi-modal SNP data 

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/6_snp_visualization/1_ArchRBrowserWithExternalData.R)

## To prepare for multi-modal visualization 

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/6_snp_visualization/5_plot_atac_connections.R)

## For visualization of snATAC-seq tracks and MDD SNPs

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/6_snp_visualization/6_visualize_snps.R)


## References
Granja, J.M., Corces, M.R., Pierce, S.E. et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nat Genet 53, 403–411 (2021). https://doi.org/10.1038/s41588-021-00790-6
