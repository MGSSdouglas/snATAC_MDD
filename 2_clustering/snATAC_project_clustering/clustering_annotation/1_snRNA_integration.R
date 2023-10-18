library(Seurat)
library(ArchR)
library(dplyr)
library(tidyr)
addArchRGenome("hg38")
addArchRThreads(threads = 16) 
set.seed(11)

#load snATAC-seq data
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded", force = FALSE, showLogo = FALSE)
proj

#load snRNA-seq data
seRNA <- readRDS("/home/anjali5/projects/def-cnagy/1_enhanced_harmonized_object.Rds")

#remove doublets and mixed cells from snRNA-seq
seRNA@meta.data['Cluster_combined'] <- plyr::mapvalues(seRNA@meta.data$Cluster, from = levels(seRNA@meta.data$Cluster), to = c("InPV", "InPV", "ExN", "ExN_L56", "ExN_L56", "ExN_L56",  "ExN_L56", "ExN_L56",  "ExN", "ExN","ExN_L24","ExN_L46", "ExN","ExN_L23" , "ExN_L23", "ExN_L24", "ExN_L46",  "Ast1", "Ast2", "InN_ADARB2", "InN_ADARB2", "InN_VIP", "InN_VIP", "ExN_L56", "ExN_L56",  "InN_SST", "InN_SST", "InN_LAMP5", "Mix", "OPC3", "Oli3","Oli1", "Oli2", "OPC1", "OPC2", "ExN", "ExN_L35", "ExN", "Mix", "Endo", "Micro"))
seRNA <- subset(seRNA, subset=DF.classifications!="Doublet")
seRNA
seRNA <- subset(seRNA, subset=Cluster_combined!="Mix")
seRNA
#Find variable features in snRNA-seq
seRNA <- FindVariableFeatures(seRNA, selection.method = "vst", nfeatures = 3000)

##Perform unconstrained integration with ArchR
proj <- addGeneIntegrationMatrix(ArchRProj = proj,useMatrix = "GeneScoreMatrix",matrixName ="GeneIntegrationMatrix",reducedDims = "Harmony",seRNA = seRNA,addToArrow = FALSE, groupRNA ='Broad', nameCell ="predictedCellBroad_Un",nameGroup ="predictedBroad_Un",nameScore = "predictedScoreBroad_Un", sampleCellsATAC = 20000, sampleCellsRNA = 20000, dimsToUse =1:20, nGenes = 2000)

##Find ArchR Clusters mapping (based on the abundance of labels predicted) to cell-type labels from snRNA
cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedBroad_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))

clustEx <- rownames(cM)[grep("ExN", preClust)]
clustIn <- rownames(cM)[grep("InN", preClust)]
clustOli <- rownames(cM)[grep("Oli", preClust)]
clustAst <- rownames(cM)[grep("Ast", preClust)]
clustOPC <- rownames(cM)[grep("OPC", preClust)]
clustEnd <- rownames(cM)[grep("End", preClust)]
clustMic <- rownames(cM)[grep("Mic", preClust)]

#Rename snATAC-seq clusters based on snRNA-seq labels mapped
cM <- confusionMatrix(proj$Clusters, proj$predictedBroad_Un)
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
proj$ClustersMapped <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = labelOld)
#saveArchRProject(proj)
#rm(cM)
#rm(preClust)

#Get snRNA-seq barcodes for each cell-type
rnaEx <- colnames(seRNA)[grep("ExN", seRNA$Cluster)]
rnaIn <- colnames(seRNA)[grep("InN", seRNA$Cluster)]
rnaOli <- colnames(seRNA)[grep("Oli", seRNA$Cluster)]
rnaOPC <- colnames(seRNA)[grep("OPC", seRNA$Cluster)]
rnaAst <- colnames(seRNA)[grep("Ast", seRNA$Cluster)]
rnaEnd <- colnames(seRNA)[grep("End", seRNA$Cluster)]
rnaMic <- colnames(seRNA)[grep("Mic", seRNA$Cluster)]

#create grouped list for constrained integration b/w snATAC-seq and snRNA-seq
groupList <- SimpleList(Ex = SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustEx],RNA = rnaEx),In = SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustIn],RNA = rnaIn), OPC=SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustOPC], RNA = rnaOPC), Oli=SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustOli], RNA = rnaOli), Ast=SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustAst], RNA = rnaAst), End=SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustEnd], RNA = rnaEnd), Mic=SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustMic], RNA = rnaMic))

#constrained integration to impute snRNA-seq values
#proj <- addGeneIntegrationMatrix(ArchRProj = proj,useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix",reducedDims = "Harmony",seRNA = rna, addToArrow = FALSE, groupRNA = 'Cluster',nameCell = "predictedCellGroup_Cn",nameGroup = "predictedGroup_Cn",nameScore = "predictedScoreGroup_Cn",sampleCellsATAC = 20000, sampleCellsRNA = 20000, dimsToUse = 1:20, nGenes = 2000)
proj <- addGeneIntegrationMatrix(ArchRProj = proj,useMatrix = "GeneScoreMatrix",matrixName ="GeneIntegrationMatrix",reducedDims = "Harmony", seRNA = seRNA,addToArrow = TRUE, groupRNA ='Broad',nameCell ="predictedCellBroad_Cn",nameGroup = "predictedBroad_Cn",nameScore = "predictedScoreBroad_Cn", sampleCellsATAC = 20000, sampleCellsRNA = 20000, dimsToUse = 1:20, nGenes = 2000, groupList=groupList,force=TRUE)
#print("done")
saveArchRProject(ArchRProj = proj, load=TRUE)
#print("saved")

cM <- as.matrix(confusionMatrix(proj$ClustersMapped, proj$predictedBroad_Un))
cM
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))

df <- proj@cellColData
dplyr::count(data.frame(df), predictedScoreBroad_Un < 0.5)

#plot histogram of unconstrained predicted lables
pdf("/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/histogram_predictedScore_broad.pdf")
hist(proj$predictedScoreBroad_Un)
abline(v = median(proj$predictedScoreBroad_Un), col="red")
dev.off()

#saveArchRProject(ArchRProj = proj, load=TRUE)
sessionInfo()
