library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 32) 
set.seed(12)


proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_FinePeaks_Split_Filtered")
#marker peaks for each cluster
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "SubClusters", #ClustersMapped for broad
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & Log2FC > 0.5")
lapply(markerList,function(x){print(dim(x))})
#write to output
saveRDS(markerList, "~/scratch/cluster_marker_peaks.rds")
#plot
pdf("/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_Peaks/Plots/Peak=Heatmap")
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR < 0.05 & Log2FC > 0.5",
  transpose = TRUE,returnMatrix = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
