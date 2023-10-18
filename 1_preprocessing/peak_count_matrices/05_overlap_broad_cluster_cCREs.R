library(Seurat)
library(Signac)
library(IRanges)

DefaultAssay(obj1) <- "peaks"
DefaultAssay(obj2) <- "peaks"

obj1 <- readRDS("ArchR_cluster_by_peaks_updated.rds") #cluster peak count matrix
obj2 <- readRDS("snATAC_broadpeaks_seuObj.rds")  # broad peak count matrix

peaks1 <- Signac::StringToGRanges(rownames(obj1))
peaks2 <- Signac::StringToGRanges(rownames(obj2))

#test overlap of broad and cluster level peaks 

length(subsetByOverlaps(peaks1, peaks2)) / length(peaks1)
#[1] 0.6388202
length(subsetByOverlaps(peaks2, peaks1)) / length(peaks2)
#[1] 0.8550953

length(subsetByOverlaps(peaks1, peaks2,maxgap=5000)) / length(peaks1)
#[1] 0.9525883
