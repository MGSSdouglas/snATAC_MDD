#load required paclages 
library(ArchR)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(cowplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
#library(spatialLIBD)
#library(GenomeInfoDb)
#library(GenomicRanges)

#######load and renormalize atac dataset

Dir_save="/home/anjali5/scratch/graham_merge_analysis/merged_Rdata/"
ArchR <- readRDS("/home/anjali5/scratch/cluster_profile/Mic_snRNA_GeneIntegrationMatrix_ArchR.rds")
ArchR
Idents(ArchR) <- "SubClusters"
ArchR <- subset(ArchR, idents=c("Mic1","Mic2"))
table(ArchR$SubClusters)
ArchR <- SCTransform(ArchR, assay = "originalexp") %>% RunPCA(verbose = FALSE, npcs = 50)
ArchR@active.assay <- "SCT"

#load Spatial data
spatial_male <- readRDS("/home/anjali5/projects/def-cnagy/anjali5/Male_merged_sections_sptial_LIBD.Rds")
DefaultAssay(spatial_male) <- "Spatial"
##female data
spatial_female <- readRDS("/home/anjali5/projects/def-cnagy/anjali5/Female_merged_sections_sptial_LIBD.Rds")
DefaultAssay(spatial_female) <- "Spatial"
#merge male and female data
spatial <- merge(x= spatial_male,y=spatial_female)
spatial <- SCTransform(spatial, assay = "Spatial", verbose = TRUE) %>% RunPCA(npcs = 50,verbose = FALSE) %>% RunUMAP(dims = 1:30)         
#find transfer anchors and make predictions 
anchors <- FindTransferAnchors(reference = ArchR, query = spatial, normalization.method = "SCT",npcs = 30)
predictions.assay <- TransferData(anchorset = anchors, refdata = ArchR$SubClusters, prediction.assay = TRUE, weight.reduction = spatial[["pca"]], dims = 1:30)
spatial[["predictions"]] <- predictions.assay
saveRDS(spatial,paste0("/home/anjali5/scratch/SpatialVisualization/mergedSections_Mic_GeneExp_Spatial.rds"))
print("DONE")
