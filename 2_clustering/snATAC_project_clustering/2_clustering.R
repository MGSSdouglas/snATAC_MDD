##Load libraries##########################
library(ArchR)
library(Seurat)
library(Signac)
library(S4Vectors)
library(dplyr)
library(tidyr)

proj <- loadArchRProject("/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded")

###########Clustering########################
proj <- addClusters(input = proj,reducedDims = "Harmony", method = "Seurat",name = "Clusters", resolution = 0.4, algorithm=3, dimsToUse = 1:25, seed=3, filterBias = TRUE)
proj <- addUMAP(ArchRProj = proj,reducedDims = "Harmony", name = "UMAPHarmony", metric = 'cosine',dimsToUse = 1:25, seed=41, n_neighbors=30,force=TRUE)
proj <- addTSNE(ArchRProj = proj,reducedDims = "Harmony",name = "TSNEHarmony",dimsToUse = 1:25, seed=5, perplexity=30,force=TRUE)

##save########
#proj <- proj4
#p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
#p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
#plotPDF(p1,p2, name = "UMAP-Proj-Filtered",ArchRProj =proj)

saveArchRProject(ArchRProj = proj,load=TRUE)
print("saved ArchR project")

#####Cluster QC#####
#proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded")
#proj

##Print cells in every Clusters
df <- proj@cellColData
with(df, table(Sample,Clusters))
table(proj$Clusters)


##Sub-clustering excitatory neurons #"C7", "C8", "C9"##
proj2 <- proj[proj$Clusters %in% c("C7", "C8", "C9"),]
proj2 <- subsetArchRProject(proj2,cells = getCellNames(proj2),dropCells = TRUE,force=TRUE)
proj2 <- addIterativeLSI(proj2, force=TRUE, dimsToUse=1:20,varFeatures = 15000, iterations=2, clusterParams = list(resolution = c(0.2)))
proj2 <- addHarmony(ArchRProj = proj2, reducedDims = "IterativeLSI",name = "Harmony",groupBy = "batch", max.iter.cluster = 20, force=TRUE)
proj2 <- addClusters(input = proj2,reducedDims = "Harmony", method = "Seurat",name = "Clusters",resolution = 0.6, knnAssign = 30, nOutlier = 50,filterBias = TRUE,algorithm=3,dimsToUse = 1:20, seed=3, force=TRUE)
proj2 <- addUMAP(ArchRProj = proj2,reducedDims = "Harmony", name = "UMAPHarmony", metric = 'cosine', dimsToUse=1:20,force=TRUE)
proj2 <- addTSNE(ArchRProj = proj2,reducedDims = "Harmony",name = "TSNEHarmony",dimsToUse = 1:20,force=TRUE)
saveArchRProject(proj2,outputDirectory  = "/home/anjali5/scratch/ArchR_outputs/ArchR_Neurons",dropCells = TRUE,load = FALSE)

##Sub-clustering Oligodendrocytes #"C1","C2","C3"##
proj_sub <- proj[proj$Clusters %in% c("C1","C2","C3"),]
proj_sub <- addIterativeLSI(proj_sub, name = "IterativeLSI", iterations = 2,varFeatures = 15000,clusterParams = list(resolution = c(0.2)),dimsToUse = 1:20,force=TRUE, seed=249)
proj_sub <- addHarmony(ArchRProj = proj_sub, reducedDims = "IterativeLSI",name = "Harmony",groupBy = "batch",max.iter.cluster = 20,dimsToUse = 1:20,force=TRUE)
proj_sub <- addClusters(input = proj_sub,reducedDims = "Harmony",method = "Seurat",resolution = 0.6, algorithm=3,knnAssign = 30, nOutlier = 50,dimsToUse = 1:20, filterBias = TRUE, force=TRUE, name="Clusters",seed=12)
proj_sub <- addUMAP(ArchRProj = proj_sub,reducedDims = "Harmony", name = "UMAPHarmony", metric = 'cosine', dimsToUse=1:20,,force=TRUE)
proj_sub <- addTSNE(ArchRProj = proj_sub,reducedDims = "Harmony",name = "TSNEHarmony",dimsToUse = 1:20,force=TRUE)
saveArchRProject(ArchRProj = proj_sub,outputDirectory ="/home/anjali5/scratch/ArchR_outputs/ArchR_Oli",,dropCells = TRUE,load = FALSE)

print("done")
sessionInfo()

