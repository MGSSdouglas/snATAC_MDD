library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 32) 
set.seed(12)

#Add Broad & Fine peaks with reproducibility=2
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded", force = FALSE, showLogo = FALSE)
proj2 <- proj
#remove low quality clusters from peak calling
proj2 <- proj2[proj2$SubClusters %ni% c("InN1", "InN2"),]
proj2
rm(proj)
#create four pseudobulk replicates of cells for peak-calling #each representing both sexes and conditions
df <- proj2@cellColData
df$sample_split <- ifelse(df$Sample %in% c("atac_A1B1", "atac_A2B2", "atac_A3B3", "atac_A4B4", "atac_A25B25", "atac_A26B26", "atac_A27B27", "atac_A28B28", "atac_A29B29"),"1",ifelse(df$Sample %in% c("atac_A9B9", "atac_A10B10","atac_A30B30","atac_A31B31","atac_A32B32", "atac_A33B33","atac_A34B34", "atac_A15B15", "atac_A16B16", "atac_A17B17"), "2", ifelse(df$Sample %in% c("atac_A35B35", "atac_A36B36", "atac_A37B73","atac_A38B38","atac_A39B39","atac_A22B22","atac_A11B11","atac_A12B12","atac_A13B13", "atac_A20B20","atac_A41B41", "atac_A42B42"), "3", "4")))
proj2@cellColData <- df
proj2
saveArchRProject(ArchRProj = proj2, outputDirectory = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered", load=FALSE)
saveArchRProject(ArchRProj = proj2, outputDirectory = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered", load=FALSE)
rm(proj2)
rm(list=ls())
#call peaks with MACS2 in broad cell-types
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered")
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersMapped",sampleLabels = "sample_split",useLabels = TRUE,minCells =683, maxCells =20588, minReplicates = 4,maxReplicates = 4,force = TRUE)
proj <- addReproduciblePeakSet(ArchRProj = proj,groupBy = "ClustersMapped",pathToMacs2="/home/anjali5/.local/bin/macs2", reproducibility="2", minCells = 683, excludeChr = "chrM",force=TRUE)
proj <- addPeakMatrix(proj)
saveArchRProject(ArchRProj = proj,load=FALSE)
rm(proj)
rm(list=ls())
#call peaks with MACS2 in clusters
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "SubClusters",sampleLabels = "sample_split",useLabels = TRUE,minCells =112, maxCells =6908, minReplicates = 4,maxReplicates = 4,force = TRUE)
proj <- addReproduciblePeakSet(ArchRProj = proj,groupBy = "SubClusters",pathToMacs2="/home/anjali5/.local/bin/macs2", reproducibility="2", minCells =112, excludeChr = "chrM", force=TRUE)
proj <- addPeakMatrix(proj)
saveArchRProject(ArchRProj = proj)

print("added peaks and saved")
#rm(list=ls())


