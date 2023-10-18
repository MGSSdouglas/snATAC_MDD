library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 32) 
set.seed(12)

proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_FinePeaks_Split_Filtered")

#Enrichment of FANS bulk ATAC Hauberg et.al Peaks 
#https://www.nature.com/articles/s41467-020-19319-2
FansPeaks <- c(GABA="~/scratch/GABA.bed",GLUT="~/scratch/GLUT.bed",MGAS="~/scratch/MGAS.bed",OLIG="~/scratch/OLIG.bed")
proj <- addPeakAnnotations(ArchRProj = proj, regions = FansPeaks, name = "FANSPeaks")
enrichRegions <- peakAnnoEnrichment(seMarker = markersPeaks,ArchRProj = proj,peakAnnotation = "FANSPeaks",cutOff = "FDR < 0.05 & Log2FC > 0.5")
heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 1, transpose = TRUE,returnMatrix = TRUE)

#Plot
print("Plotting..")
pdf("~/scratch/FullardPeaks_In_ClusterMarkers.pdf",height=10,width=8)
ComplexHeatmap::Heatmap(matrix = heatmapRegions,col = paletteContinuous(set = "comet", n = 100),column_names_gp = grid::gpar(fontsize = 15),row_names_gp = grid::gpar(fontsize = 15),cluster_rows = TRUE,cluster_row_slices = TRUE,column_title = "ATAC_FANS_PEAKS",row_title = "snATAC",clustering_distance_rows ="binary",clustering_distance_columns = "binary",heatmap_legend_param = list(title = "Enrichment"))
dev.off()

