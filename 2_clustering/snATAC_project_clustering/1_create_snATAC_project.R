##Load libraries##########################
library(ArchR)
library(Seurat)
library(Signac)
library(S4Vectors)
library(dplyr)
library(tidyr)

##Add genome###############################
addArchRGenome("hg38")
addArchRThreads(threads = 32)
set.seed(1234)

##Create ArchR arrow files with demultiplexed subjects######################
inputFiles = c("atac_A1B1"="/home/anjali5/scratch/graham_atac_analysis/atac_A1B1/fragments.tsv.gz", "atac_A2B2"="/home/anjali5/scratch/graham_atac_analysis/atac_A2B2/fragments.tsv.gz","atac_A3B3"="/home/anjali5/scratch/graham_atac_analysis/atac_A3B3/fragments.tsv.gz","atac_A4B4"="/home/anjali5/scratch/graham_atac_analysis/atac_A4B4/fragments.tsv.gz","atac_A5B5"="/home/anjali5/scratch/graham_atac_analysis/atac_A5B5/fragments.tsv.gz","atac_A6B6"="/home/anjali5/scratch/graham_atac_analysis/atac_A6B6/fragments.tsv.gz","atac_A7B7"="/home/anjali5/scratch/graham_atac_analysis/atac_A7B7/fragments.tsv.gz","atac_A8B8"="/home/anjali5/scratch/graham_atac_analysis/atac_A8B8/fragments.tsv.gz","atac_A9B9"="/home/anjali5/scratch/graham_atac_analysis/atac_A9B9/fragments.tsv.gz","atac_A10B10"="/home/anjali5/scratch/graham_atac_analysis/atac_A10B10/fragments.tsv.gz","atac_A11B11"="/home/anjali5/scratch/graham_atac_analysis/atac_A11B11/fragments.tsv.gz","atac_A12B12"="/home/anjali5/scratch/graham_atac_analysis/atac_A12B12/fragments.tsv.gz","atac_A13B13"="/home/anjali5/scratch/graham_atac_analysis/atac_A13B13/fragments.tsv.gz","atac_A14B14"="/home/anjali5/scratch/graham_atac_analysis/atac_A14B14/fragments.tsv.gz","atac_A15B15"="/home/anjali5/scratch/graham_atac_analysis/atac_A15B15/fragments.tsv.gz","atac_A16B16"="/home/anjali5/scratch/graham_atac_analysis/atac_A16B16/fragments.tsv.gz","atac_A17B17"="/home/anjali5/scratch/graham_atac_analysis/atac_A17B17/fragments.tsv.gz","atac_A18B18"="/home/anjali5/scratch/graham_atac_analysis/atac_A18B18/fragments.tsv.gz","atac_A19B19"="/home/anjali5/scratch/graham_atac_analysis/atac_A19B19/fragments.tsv.gz","atac_A20B20"="/home/anjali5/scratch/graham_atac_analysis/atac_A20B20/fragments.tsv.gz","atac_A21B21"="/home/anjali5/scratch/graham_atac_analysis/atac_A21B21/fragments.tsv.gz","atac_A22B22"="/home/anjali5/scratch/graham_atac_analysis/atac_A22B22/fragments.tsv.gz","atac_A23B23"="/home/anjali5/scratch/graham_atac_analysis/atac_A23B23/fragments.tsv.gz","atac_A24B24"="/home/anjali5/scratch/graham_atac_analysis/atac_A24B24/fragments.tsv.gz","atac_A25B25"="/home/anjali5/scratch/graham_atac_analysis/atac_A25B25/fragments.tsv.gz","atac_A26B26"="/home/anjali5/scratch/graham_atac_analysis/atac_A26B26/fragments.tsv.gz","atac_A27B27"="/home/anjali5/scratch/graham_atac_analysis/atac_A27B27/fragments.tsv.gz","atac_A28B28"="/home/anjali5/scratch/graham_atac_analysis/atac_A28B28/fragments.tsv.gz","atac_A29B29"="/home/anjali5/scratch/graham_atac_analysis/atac_A29B29/fragments.tsv.gz","atac_A31B31"="/home/anjali5/scratch/graham_atac_analysis/atac_A31B31/fragments.tsv.gz","atac_A32B32"="/home/anjali5/scratch/graham_atac_analysis/atac_A32B32/fragments.tsv.gz","atac_A33B33"="/home/anjali5/scratch/graham_atac_analysis/atac_A33B33/fragments.tsv.gz","atac_A34B34"="/home/anjali5/scratch/graham_atac_analysis/atac_A34B34/fragments.tsv.gz","atac_A35B35"="/home/anjali5/scratch/graham_atac_analysis/atac_A35B35/fragments.tsv.gz","atac_A37B37"="/home/anjali5/scratch/graham_atac_analysis/atac_A37B37/fragments.tsv.gz","atac_A38B38"="/home/anjali5/scratch/graham_atac_analysis/atac_A38B38/fragments.tsv.gz","atac_A39B39"="/home/anjali5/scratch/graham_atac_analysis/atac_A39B39/fragments.tsv.gz","atac_A40B40"="/home/anjali5/scratch/graham_atac_analysis/atac_A40B40/fragments.tsv.gz","atac_A41B41"="/home/anjali5/scratch/graham_atac_analysis/atac_A41B41/fragments.tsv.gz", "atac_A42B42"="/home/anjali5/scratch/graham_atac_analysis/atac_A42B42/fragments.tsv.gz","atac_A30B30"="/home/anjali5/scratch/graham_atac_analysis/atac_A30B30/fragments.tsv.gz","atac_A36B36"="/home/anjali5/scratch/graham_atac_analysis/atac_A36B36/fragments.tsv.gz")
ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = names(inputFiles),minTSS = 3,minFrags = 1000,addTileMat = TRUE,addGeneScoreMat = TRUE)
ArrowFiles

##Add ArchR Doublet scores#################
doubScores <- addDoubletScores(input = ArrowFiles, k = 10,knnMethod = "UMAP")

###Create ArchR Project#####################
proj <- ArchRProject(ArrowFiles = ArrowFiles,  outputDirectory = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded",copyArrows = FALSE)
proj

##Add batch information######################
l1=c("atac_A5B5","atac_A7B7", "atac_A17B17", "atac_A25B25")
l2=c("atac_A29B29","atac_A34B34","atac_A35B35","atac_A39B39")
l3=c("atac_A1B1","atac_A6B6","atac_A31B31","atac_A40B40")
l4=c("atac_A2B2","atac_A8B8","atac_A24B24","atac_A37B37")
l5=c("atac_A3B3","atac_A18B18","atac_A19B19","atac_A33B33")
l6=c("atac_A4B4","atac_A12B12","atac_A22B22", "atac_A30B30")
l8=c("atac_A9B9","atac_A10B10","atac_A26B26","atac_A41B41")
l9=c("atac_A11B11","atac_A13B13","atac_A15B15","atac_A21B21","atac_A27B27")
l10=c("atac_A16B16","atac_A28B28", "atac_A32B32","atac_A38B38")
l11=c("atac_A14B14","atac_A20B20","atac_A23B23", "atac_A36B36")

#Adding batch to metadata####################
df <- proj@cellColData
df$batch <- ifelse(df$Sample %in% l1,'1', ifelse(df$Sample %in% l2, '2', ifelse(df$Sample %in% l3, '3', ifelse(df$Sample %in% l4, '4', ifelse(df$Sample %in% l5, '5', ifelse(df$Sample %in% l6, '6',ifelse(df$Sample %in% l8, '8', ifelse(df$Sample %in% l9, '9',ifelse(df$Sample %in% l10, '10',ifelse(df$Sample %in% l11,'11','12'))))))))))
with(df, table(batch,Sample))
proj@cellColData <- df

##Add more metadata to ArchR#############
Dir="/home/anjali5/scratch/graham_atac_analysis/"
df <- proj@cellColData
dim(df)
foo <- data.frame(do.call('rbind', strsplit(as.character(rownames(df)),'#',fixed=TRUE)))
df$barcode <- foo$X2
df$row <- rownames(df)
meta <- data.frame()
add_metadata <- function(meta, val){
  print(val)
  tmp <- readRDS(paste0(Dir,"atac_",val, "/atac_sex.rds"))
  tmp
  md <- tmp@meta.data
  md$Sample <- paste0("atac_", val)
  meta <- rbind(meta,md)
  dim(meta)
  return(meta)
}

list=c("A1B1", "A2B2", "A3B3", "A4B4", "A5B5", "A6B6", "A7B7", "A8B8", "A9B9", "A10B10", "A11B11", "A12B12", "A13B13", "A14B14", "A15B15", "A16B16", "A17B17", "A18B18", "A19B19", "A20B20", "A21B21", "A22B22", "A23B23", "A24B24", "A25B25", "A26B26", "A27B27", "A28B28", "A29B29", "A31B31", "A32B32", "A33B33", "A34B34", "A35B35",  "A37B37", "A38B38", "A39B39", "A40B40", "A41B41", "A42B42", "A30B30", "A36B36")
for(i in list){
  print(i)
  meta <- add_metadata(meta,i)
}

dim(meta)
dim(df)
#Merge metadata with ArchR project######################
df_meta <- merge(meta,df,by=c("Sample","barcode"))
dim(df_meta)
head(df_meta)
rownames(df_meta) <- df_meta$row
typeof(df_meta)
write.csv(df_meta, "/home/anjali5/scratch/df_meta", quote=FALSE, row.names=FALSE)
#df <- as(df, "DataFrame")
#typeof(df)
proj@cellColData <- df_meta

####sanity check#######################
#df <- proj@cellColData
#df <- data.frame(df)
#dplyr::count(df, Sample=="atac_A10B10" & sex=="male")
#dplyr::count(df, Sample=="atac_A25B25" & cellsnp_doublet=="1")
#rm(df)

####Add subject demographics and phenotype information##########
my_data <- read.csv("~/scratch/sample_phenotype.xlsx")
row= c("42", "A42+B42", "M", "atac_A42B42", "Case")
new_data <- rbind(my_data,row)
case_list <- dplyr::filter(new_data, condition=="Case")
case_list <- case_list$Sample
df <- proj@cellColData
df[['condition']] <- rep("unknown", nrow(df))
df$condition <- ifelse(df$Sample %in% case_list, "Case", "Control")
proj@cellColData <- df

######Filter Low-quality cells based on TSSEnrichment & Mitochondrial contamination##########
idxSample <- BiocGenerics::which(proj$TSSEnrichment>3.5)
cellsSample <- proj$cellNames[idxSample]
proj <- proj[cellsSample, ]
proj
##
df <- proj@cellColData
df$mito_percent <- df$mitochondrial / df$total
proj@cellColData <- df
idxSample <- BiocGenerics::which(proj$mito_percent<0.1)
cellsSample <- proj$cellNames[idxSample]
proj <- proj[cellsSample, ]
proj

###Remove cells with unassigned sex###########################
idxSample <- BiocGenerics::which(proj$sex %in% c("male","female"))
cellsSample <- proj$cellNames[idxSample]
proj2 <- proj[cellsSample, ]
proj2
df <- proj2@cellColData
head(df)
rm(df)
#rm(proj)

#####Add sex information to subject name#####################
df <- proj2@cellColData
df$Subject <-  paste0(df$Sample, df$sex)
proj2@cellColData <- df

###Remove Doublets####################################
# Vireo mixed genotypes & 10x cellranger doublets########
idxSample <- BiocGenerics::which(proj2$cellsnp_doublet=="unknown" & proj2$excluded_reason=="0" | proj2$cellsnp_doublet=="unknown" & proj2$excluded_reason=="2")
cellsSample <- proj2$cellNames[idxSample]
proj3 <- proj2[cellsSample, ]
proj3
rm(proj2)

####Amulet & ArchR doublets####
#proj3 <- filterDoublets(proj3, filterRatio=1)
#proj3
#idxSample <- BiocGenerics::which(proj3$amulet!="1")
#cellsSample <- proj3$cellNames[idxSample]
#proj4 <- proj3[cellsSample, ]
#proj4
#rm(proj3)
####
proj <- proj3
rm(proj3)

######Dimensionality reduction with Iterative LSI#############################
proj <- addIterativeLSI(ArchRProj = proj, seed=1,dimsToUse = 1:30, varFeatures = 25000, clusterParams = list(resolution = c(2)))
#proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded")

#####Batch-correction##############
proj <- addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI",name = "Harmony",groupBy = "batch", dimsToUse = 1:30, max.iter.cluster = 20)

####Get ArchR GeneScore based markers#########
markersGS <- getMarkerFeatures(ArchRProj = proj,useMatrix = "GeneScoreMatrix",groupBy ="Clusters",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR < 0.05 & Log2FC > 0.5")

##Define some marker genes######################
#markerGenes  <- c("RBFOX3", "SNAP25", "SATB2","SLC17A7", "SLC17A6", "NEUROD6", "GAD1", "GAD2","SLC32A1","SLC1A2","GFAP","GLUL","AQP4","SOX9","ALDH1A1","VIM","CSF1R","CX3CR1","CD14","MRC1","TMEM119","CLDN5", "VTN", "FLT1", "SRGN", "P2RY12", "MBP","MOG","MAG","OPALIN","OLIG1","OLIG2","CSPG4","ASCL1", "HAS2","PCDH15", "LHX6","PVALB", "SST", "ADARB2","VIP", "LAMP5", "SNCG", "CUX2","RORB","NTNG2","BCL11B","THEMIS", "FEZF2")

###Impute GeneScore with MAGIC##############
#proj <- addImputeWeights(proj)
#saveArchRProject(ArchRProj = proj,load=TRUE)

##Plot pseudo-bulk gene covergae###############
p1 <- plotBrowserTrack(ArchRProj = proj,groupBy = "Clusters",geneSymbol = markerGenes, upstream = 5000, downstream = 5000)
plotPDF(plotList = p1, name = "Plot-Tracks-BroadClusters-Marker-Genes.pdf",ArchRProj = proj,addDOC = FALSE, width = 8, height = 10)
#saveArchRProject(ArchRProj = proj,load=TRUE)







