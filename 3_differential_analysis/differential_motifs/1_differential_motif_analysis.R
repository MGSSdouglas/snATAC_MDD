library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRGenome("hg38")
addArchRThreads(threads = 16)
set.seed(12)


proj <- loadArchRProject("~/scratch/ArchR_backup/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")

#Add Motif
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif",force=TRUE)
saveArchRProject(proj)

#Add ChromVar Motif matrix
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(proj)
print("ChromVar done")

###Marker TF motifs#################
total_df <- data.frame()
markersGS <- getMarkerFeatures(ArchRProj = proj,testMethod = "wilcoxon",groupBy="SubClusters",binarize = FALSE,useMatrix = "MotifMatrix",useSeqnames="z",bias = c("log10(nFrags)","TSSEnrichment"),maxCells=1000)
markerList <- getMarkers(markersGS, cutOff = "FDR < 0.05")
for(name in names(markerList)){
  tmp <- markerList[[name]] %>% data.frame()
  tmp$cluster_id <- name
  total_df <- rbind(total_df,tmp)
}
write.csv(total_df, "~/projects/def-gturecki/anjali5/subcluster_marker_tf_result.csv")

###Case vs control TF motifs#########

proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")
proj$PMI <- as.numeric(proj$PMI)
proj$Age <- as.numeric(proj$Age)

total_df <- data.frame()

for(cluster in unique(proj$SubClusters)){ 
  print(cluster)
  proj2 <- proj[proj$SubClusters %in% cluster,]
  diffMotif <- getMarkerFeatures(ArchRProj = proj2,testMethod = "wilcoxon",useGroups = "Case", bgdGroups = "Control",binarize = FALSE,useMatrix = "MotifMatrix",groupBy = "condition",useSeqnames="z",maxCells = 1000,bias = c("log10(nFrags)", "PMI", "Age"),threads = 16)
  markerList <- getMarkers(diffMotif, cutOff = "FDR < 0.05")
  markerList$Case %>% data.frame() -> tmp_df
  tmp_df$cluster_id <- rep(cluster,nrow(tmp_df))
  total_df <- rbind(total_df, tmp_df)
  print(dim(total_df))
}
write.csv(total_df, "~/projects/def-gturecki/anjali5/diff_tfs_chromvar.csv")

##Plot TF footprints####

#first load edited TF footprint code
#https://github.com/GreenleafLab/ArchR/issues/493
source("Get_ArchR_Footprints.R")

motifPositions <- getPositions(proj, name="Motif")

motifs <- c("SOX9", "JUN","SPI1","LHX2","BHLHA15","JUNB","FOSL2","DLX5", "POU2F2", "POU3F4","NEUROD1","SPIB", "NFIA")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClustersMapped",force=T)
seFoot <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[markerMotifs],
  groupBy = "ClustersMapped"
)
plotFootprintsManual(   #the result will store in ArchR Project folder under "Plots"
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-Broad-Marker-Colored", #name this unique
  addDOC = FALSE,
  smoothWindow = 5, pal = c("Ast"="#145F04","End"="#E3D2A7", "ExN"="#b36609", "InN"="#4853A2", "Mic"="#1792AA", "Oli"="#C60E9D", "OPC"="#f08080")
)


proj$SubClusters <- plyr::mapvalues(proj$SubClusters, from = unique(proj$SubClusters), to = c("ExN12","Ast1", "ExN5", "ExN9", "InVIP", "ExN10", "Ast3", "ExN2","Oli5", "ExN6", "Mix2", "InSST", "ExN3", "ExN11", "InPV", "ExN8", "ExN7", "Oli6", "Oli7", "End1", "ExN1", "InLAMP5", "Ast4", "OPC2", "Oli4", "Mic2", "InN3", "Oli2","Ast2", "ExN3", "Oli1", "Mic1", "OPC4", "OPC1", "OPC3", "Mix1", "Oli3", "End2"))

motifs <- c("JUN","FOSL2","EGR1","NEUROD1", "POU2F2", "POU3F4", "SOX9", "SPI1", "SPIB", "NFIA", "NFIB", "IRF1","NR4A2")

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "SubClusters",force=T)
seFoot <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[markerMotifs],
  groupBy = "SubClusters"
)

col= c("Ast1"="#23B205","Ast2"="#74C69D","Ast3"="#2DF304","Ast4"="#145F04","End1"="#E3D2A7","End2"="#F7E3B4","ExN1"="#E
7B321","ExN11"="#A6A3A1","ExN10"="#CCC7C4","ExN7"="#D0C75C","ExN3"="#FAB866","ExN2"="#FAE605","ExN12"="#736F6C", "ExN8"="#F7F7B5","ExN4"="#F29F3A","ExN9"="#E7D973","ExN5"="#D18221","ExN6"="#B36609","InLAMP5"="#53B1F3","InPV"="#4853A2","InSST
"="#6E7EBD", "InVIP"="#81B2F1","InN3"= "#C0D9FC", "Mic1" = "#1792AA", "Mic2" = "#83D6E6", "Mix1" = "#D2D6CB", "Mix2" = "#D7C7D5", "Oli1" = "#C60E9D","Oli2"="#B199C7","Oli3"= "#C3B2D2","Oli4"= "#7F56A4","Oli5"= "#E04EC0","Oli6"= "#E677CD","Oli7"= "#F4BBE7", "OPC1"= "#F4978E", "OPC2"="#FFB4A2", "OPC3"= "#F8651B", "OPC4"= "#F08080")
plotFootprintsManual(   #the result will store in ArchR Project folder under "Plots"
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-Sub-Marker-Colored", #name this unique
  addDOC = FALSE,
  smoothWindow = 5,pal=col)
}

#case-control motifs in ExN1 and Mic2
#ExN1
proj2 <- proj[proj$SubClusters %in% "ExN1",]
motifPositions <- getPositions(proj, name="Motif")
motifs <- c("BACH1", "BACH2", "FOSL2", "CTCF", "EGR1", "NEUROD1", "ASCL1", "RFX")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

#10rep
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "condition",maxReplicate=10,force=T)
seFoot <- getFootprints( ArchRProj = proj2,positions = motifPositions[markerMotifs], groupBy ="condition")
plotFootprintsManual(seFoot = seFoot,ArchRProj = proj,plotName = "Cisbp-Updated-10sub-ExN1-TFs-Footprints-Subtract-Bias")

#Mic2
proj2 <- proj[proj$SubClusters %in% "Mic2",]
motifs <- c("SP","KLF","FLI1","IRF","ETV","CRE","JUND","ATF","CREM","RUNX","SMAD")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

#10rep
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "condition",maxReplicate=10,force=T)
seFoot <- getFootprints( ArchRProj = proj2,positions = motifPositions[markerMotifs], groupBy ="condition")
plotFootprintsManual(seFoot = seFoot,ArchRProj = proj,plotName = "Cisbp-Updated-10sub-Mic2-TFs-Footprints-Subtract-Bias")


