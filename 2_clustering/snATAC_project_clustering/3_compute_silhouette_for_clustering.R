##Load libraries
library(ArchR)
library(Seurat)
library(Signac)
library(readr)
library(dplyr)
library(ggplot2)

addArchRGenome("hg38")
addArchRThreads(threads = 32)
set.seed(2) 

####Silhoutte##########
CalculateSilhouette<- function(proj,clust){
  library(cluster)
  library(readr)
  rd <- getReducedDims(proj,reducedDims = "Harmony")
  cell_distance<- dist(rd)
  print(clust)
  cell_cluster<- parse_number(getCellColData(proj, select = clust, drop = TRUE))
  silhouette_score <- cluster::silhouette(cell_cluster, cell_distance)
  proj <- addCellColData(ArchRProj = proj, data = paste0(silhouette_score[, 3]),name = paste0("silhouette_",res), cells=proj$cellNames,force=TRUE)
  #proj$paste0("silhouette_",res)  <-  silhouette_score[, 3]
  return(proj)
}
######Sub-Clustering####

proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded")
clusters_to_subcluster <- unique(proj$ClustersMapped)

for(cluster in clusters_to_subcluster) {
  resolutions <- c(0.6,0.4,0.2)
  #datalist = list()
  dataFrame_resolution = data.frame()
  proj_sub <- proj[proj$Clusters %in% cluster,]
  #proj_sub <- subsetArchRProject(proj_sub, cells = rownames(proj_sub),outputDirectory = paste0("/home/anjali5/scratch/ArchR_outputs/","ArchR_BatchTSS", cluster), force = TRUE)
  #proj_sub <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_C15_2")
  proj_sub <- addIterativeLSI(proj_sub, name = "IterativeLSI", iterations = 2,varFeatures = 15000,clusterParams = list(resolution = c(2)),dimsToUse = 1:20,force=TRUE, seed=249)
  proj_sub <- addHarmony(ArchRProj = proj_sub, reducedDims = "IterativeLSI",name = "Harmony",groupBy = "batch",max.iter.cluster = 20,dimsToUse = 1:20,force=TRUE)
  for(res in resolutions){
    #proj_sub <- addClusters(input = proj_sub,reducedDims = "Harmony",method = "Seurat",resolution = res, algorithm=3, nOutlier = 50, knnAssign = 30, dimsToUse = 1:20, filterBias = TRUE, force=TRUE, name=paste0("Clusters_",res),seed=12)
    proj_sub <- addClusters(input = proj_sub,reducedDims = "Harmony",method = "Seurat",resolution =res, algorithm=3, k.param=30, dimsToUse = 1:20, nOutlier=50, knnAssign = 30, filterBias = TRUE, force=TRUE, name=paste0("Clusters_",res),seed=12)
  }
  for(res in resolutions){
    data= data.frame()
    tmp <- paste0("Clusters_",res)
    if(data.frame(proj_sub@cellColData) %>% summarise(count = n_distinct(get(tmp))) < 2){
      print(paste0("only one cluster found at:", res))
      #print("Do not sub-cluster unless resolution less than current is promising..")
      print(paste0("not sufficient clusters to calculate sil at:", res))
      break
    }
    print(paste0("calculating silhoutee .. at resolution:",res))
    proj_sub <- CalculateSilhouette(proj_sub, tmp)
    #df <- data.frame(proj_sub@cellColData)
    tmp0 <- paste0("silhouette_",res)
    assign(paste0("Mean_",res),aggregate(as.numeric(getCellColData(proj_sub, select = tmp0, drop = TRUE)), list(getCellColData(proj_sub, select = tmp, drop = TRUE)), FUN=mean))
    assign(paste0("CellsInCluster_",res), as.data.frame(table(getCellColData(proj_sub, select = tmp))))
    print("added to dataFrame")
    data <- data.frame(resolution=res, Mean_Silh=get(paste0("Mean_",res)), cells=get(paste0("CellsInCluster_",res)), Median_Silh=median(get(paste0("Mean_",res))$x))
    print(data)
    dataFrame_resolution <- rbind(dataFrame_resolution,data)
  }
  saveArchRProject(ArchRProj = proj_sub, outputDirectory = paste0("/home/anjali5/scratch/ArchR_outputs/","ArchR_BatchTSS_", cluster))
  print(paste0("for cluster", cluster))
  print(dataFrame_resolution)
}


