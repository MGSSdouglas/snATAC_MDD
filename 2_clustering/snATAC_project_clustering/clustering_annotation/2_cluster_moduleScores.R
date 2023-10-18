##Load libraries
#library(ArchR)
library(Seurat)
library(Signac)
library(readr)
library(dplyr)
library(ggplot2)
library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 16)

##Load ArchR project
#proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters")
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded")

##Extract All ArchR Features
allFeatures <- getFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix")
#features_not_present <- features[which(features %ni% allFeatures)]

print("Adding BrainInABlender module scores")
#Add BrainInABlender module scores
library(BrainInABlender)
#This includes mouse genes 
CellTypeSpecificGenes <- dplyr::filter(CellTypeSpecificGenes_Master3, Species=="Human")
BIAB_markers <- split(CellTypeSpecificGenes$GeneSymbol_Human, CellTypeSpecificGenes$CellType_Primary)
BIAB_markers <- lapply(BIAB_markers, na.omit)
BIAB_markers <- lapply(BIAB_markers, as.character)

BIAB_markers[["Mural"]] <- NULL
BIAB_markers[["Neuron_Interneuron"]] <- NULL
BIAB_markers[["Neuron_Projection"]] <- NULL

##Identify overlapping feature-set
BIAB_markers <- lapply(BIAB_markers, function(x){
  x[which(x %in% allFeatures)]
})

##Add BIAB Module scores 
proj <- addModuleScore(proj, features = BIAB_markers, useMatrix = "GeneScoreMatrix")

##Print median module scores for very cluster
Clusters_list <- unique(proj$NewClust)
#Clusters_list <- unique(proj$SubClust)
Modules <- names(BIAB_markers)
for(cluster in Clusters_list){
  proj_tmp <- proj[proj$NewClust %in% cluster,]
  for(module in Modules){
    tmp <- paste0("Module.",module)
    x <- median(getCellColData(proj_tmp, select = tmp, drop = TRUE))
    print(paste0("For Cluster in:", cluster, "-Median score for Module-", module, " is: ", x))
  }
}

##Plot UMAP Embeddings (imputed and non-imputed)
plot_list <- list()
for(i in names(BIAB_markers)){
  plot_list[[i]] <- plotEmbedding(proj, name=paste0("Module.",i), imputeWeights = getImputeWeights(proj),embedding = "UMAPHarmony4")
}
pdf(file = paste0("/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/BIAB_Module_Score_Imputed_Magic_NewClust.pdf"),
    onefile = TRUE, height = 10, width = 12)
plot_list
dev.off()

plot_list <- list()
for(i in names(BIAB_markers)){
  plot_list[[i]] <- plotEmbedding(proj, name=paste0("Module.",i), imputeWeights = NULL,embedding = "UMAPHarmony4")
}
pdf(file = paste0("/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/BIAB_Module_Score_UnImputed_NewClust.pdf"),
    onefile = TRUE, height = 10, width = 12)
plot_list
dev.off()

print("Done..")


##Add VLMC & Pericyte Module scores 
features <- list(VLMC=c("ITIH5","ABCA9","DCN","ATP1A2","SAT1","NEAT1","SLC7A11","ABCA8","APOD","NKD1","PTN","FBLN1","LAMC3","PDGFRB","CXCL12","AKR1C2","SVIL","ABCA6","COLEC12","CYP1B1","UACA","AHNAK","TPCN1","SLC1A3","ADAM33","LAMA2","COL18A1","CYTL1","SLC6A1","EDN3","ARHGAP29","COL4A5","DAB2","BGN","LEPR","PTGDS","OGN","WDR86","PRELP","ISLR","PDGFRA","MRC2","COL15A1","SPRY4","SELENBP1","RHOJ","KIAA1755","SASH1","LHFP","RPS12P5","KANK2"), Pericyte=c("PDGFRB","IGFBP7","ATP1A2","NOTCH3","ITIH5","SLC6A12","HIGD1B","RGS5","DCN","SLC6A1","PTN","IFITM3","CALD1","GNG11","BGN","FN1","IFITM1","SLC19A1","NDUFA4L2","PLXDC1","B2M","DLC1","TXNIP","COLEC12","PRELP","ADIRF","SPARC","SLC12A7","SLC6A13","CD9","EPS8","ISYNA1","LGALS1","IFITM2","EPAS1","MYO1B","STOM","HSPB1","EMP2","HES1","CYBA","PTGDS","MCAM","MYL12A","PDE7B","CD63","NR2F2","ARHGAP29","FRZB","TNS2","GJC1"))
features <- lapply(features, function(x){ x[which(x %in% allFeatures)]})
proj<- addModuleScore(proj, features = features, useMatrix = "GeneScoreMatrix")
Clusters_list <- unique(proj$NewClust)
Modules <- names(features)
for(cluster in Clusters_list){
  proj_tmp <- proj[proj$NewClust %in% cluster,]
  for(module in Modules){
    tmp <- paste0("Module.",module)
    x <- median(getCellColData(proj_tmp, select = tmp, drop = TRUE))
    print(paste0("For Cluster in:", cluster, "-Median score for Module-", module, " is: ", x))
  }
}

saveArchRProject(proj, load=TRUE)
