closeAllConnections()
rm(list = ls())

library(ArchR)
library(Seurat)
library(SingleCellExperiment)

#h5disableFileLocking()
#HDF5_USE_FILE_LOCKING=FALSE
#RHDF5_USE_FILE_LOCKING=FALSE

#export snATAC-seq broad peak matrices to seurat object
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered")
traceback()
se <- getMatrixFromProject(proj,"PeakMatrix",threads=4)
print("matrix obtained")
se <- as(se, "SingleCellExperiment")
names(assays(se)) <- "counts"
rownames(se) <- paste0(seqnames(se), "-", start(se),"-",end(se))
seuratObj <- as.Seurat(se, data = NULL)
rD <- getReducedDims(proj)
seuratObj[["Harmony"]] <- CreateDimReducObject(embeddings = rD, key = "Harmony_", assay = DefaultAssay(seuratObj))
seuratObj
saveRDS(seuratObj, "/home/anjali5/projects/def-gturecki/anjali5/seuratObj_peakMatrix_ArchR_BroadPeaks.rds")


#export snATAC-seq cluster peak matrices to seurat object # case and control (otherwise too large!!)
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")
for(conditions in c("case","control")){
  proj2 <- proj[proj$condition %in% conditions,]
  print(table(proj2$condition)) #doublc-check
  traceback()
  se <- getMatrixFromProject(proj2,"PeakMatrix",threads=4)
  print("matrix obtained")
  se <- as(se, "SingleCellExperiment")
  names(assays(se)) <- "counts"
  rownames(se) <- paste0(seqnames(se), "-", start(se),"-",end(se))
  seuratObj <- as.Seurat(se, data = NULL)
  rD <- getReducedDims(proj)
  seuratObj[["Harmony"]] <- CreateDimReducObject(embeddings = rD, key = "Harmony_", assay = DefaultAssay(seuratObj))
  seuratObj
  saveRDS(seuratObj, paste0("/home/anjali5/projects/def-gturecki/anjali5/seuratObj_peakMatrix_ArchR_finepaks",conditions,".rds"))
}

#get gene score matrix
se <- getMatrixFromProject(proj,"GeneScoreMatrix")
se <- as(se, "SingleCellExperiment")
names(assays(se)) <- "counts"
#rownames(se) <- paste0(seqnames(se), ":", start(se),"-",end(se))
rownames(se) <- paste0(rowData(se)$name)
seuratObj <- as.Seurat(se, data = NULL)
rD <- getReducedDims(proj)
seuratObj[["Harmony"]] <- CreateDimReducObject(embeddings = rD, key = "Harmony_", assay = DefaultAssay(seuratObj))
seuratObj
saveRDS(seuratObj, "/home/anjali5/projects/def-cnagy/anjali5/seuratObj_geneMatrix_ArchR.rds")


#export microglia imputed gene-experssion matrix
proj2 <- proj[proj$ClustersMapped %in% c("Mic"),]
proj2
se <- getMatrixFromProject(proj2,useMatrix ="GeneIntegrationMatrix",threads = 4)
#rm(proj)
se <- as(se, "SingleCellExperiment")
names(assays(se)) <- "counts"
rownames(se) <- paste0(rowData(se)$name)
seuratObj <- as.Seurat(se, data = NULL)
rD <- getReducedDims(proj2)
seuratObj[["Harmony"]] <- CreateDimReducObject(embeddings = rD, key = "Harmony_", assay = DefaultAssay(seuratObj))
seuratObj
saveRDS(seuratObj, paste0("/home/anjali5/scratch/Mic_GeneExpMatrix.rds"))

