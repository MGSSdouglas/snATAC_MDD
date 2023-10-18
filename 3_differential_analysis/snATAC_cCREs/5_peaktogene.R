
library(ArchR)

#calculate peak-to-peak co-accessibility and peak-to-gene linkages for broad and cluster peaks
#broad
proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered")
#cluster
proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")

#distance parameter for linkages
dist=500000
proj <- addCoAccessibility(proj,maxDist = dist)
proj <- addPeak2GeneLinks(proj,maxDist = dist) 

#plot p2g heatmap with k-means clustering (k=10) and color by broad
p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "ClustersMapped", k = 10, varCutOffATAC = 0, varCutOffRNA = 0)
pdf("~/projects/def-cnagy/anjali5/p2gheatmapk10.pdf")
p
dev.off()

