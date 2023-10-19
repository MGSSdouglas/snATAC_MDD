## Main
This directory contains scripts for creating ArchR project and performing unsupervised clustering and annotation of >200,000 nuclei

## snATAC-seq clustering
<img width="1261" alt="Screenshot 2023-10-18 at 6 05 25 PM" src="https://github.com/MGSSdouglas/snATAC_MDD/assets/60046859/c9ccb452-c569-48a3-b372-e24951165e64">

UMAP of 201,456 Chromatin accessibility defined high-quality cells colored by 38 clusters from 7 broad cell-types. Iterative clustering within broad excitatory and inhibitory cell-types revealed 12 excitatory (ExN1-12; top), 5 inhibitory clusters (InN3/PV/SST/LAMP5/VIP; bottom). UMAP plot colored by gene-activity computed using ArchR GeneScore matrix for cortical layer (top) and interneuron lineage marker genes (bottom). Top: CUX2: Upper layer (L1-2), RORB, TOX: Middle layer (L3-4), SYNPR, NTNG2, FEZF2: Deep layer (L5-6). Bottom: LHX2, SST, PVALB: medial ganglionic eminence (MGE), ADARB2, VIP, LAMP5: caudal (CGE) lineage interneurons.

## To filter high-quality nuclei and create ArchR project

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/2_clustering/snATAC_project_clustering/1_create_snATAC_project.R)

## To choose the number of LSI dimensions
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/2_clustering/snATAC_project_clustering/clustering_qc/1_scree_plot.R)

## To perform unsupervised clustering
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/2_clustering/snATAC_project_clustering/2_clustering.R)

[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/2_clustering/snATAC_project_clustering/3_compute_silhouette_for_clustering.R)

## To perform clustering QC
[script](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/2_clustering/snATAC_project_clustering/clustering_qc)

## To perform annotation of cell-types and clusters
[script](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/2_clustering/snATAC_project_clustering/clustering_annotation)
