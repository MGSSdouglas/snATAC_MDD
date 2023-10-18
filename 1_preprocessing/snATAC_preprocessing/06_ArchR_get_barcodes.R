library(ArchR)
addArchRGenome("hg38")
addArchRThreads(threads = 32)
HDF5_USE_FILE_LOCKING=FALSE
RHDF5_USE_FILE_LOCKING=FALSE
#set.seed(1)
inputFiles = c("atac_A1B1"="/home/anjali5/scratch/graham_atac_analysis/atac_A1B1/fragments.tsv.gz", "atac_A2B2"="/home/anjali5/scratch/graham_atac_analysis/atac_A2B2/fragments.tsv.gz","atac_A3B3"="/home/anjali5/scratch/graham_atac_analysis/atac_A3B3/fragments.tsv.gz","atac_A4B4"="/home/anjali5/scratch/graham_atac_analysis/atac_A4B4/fragments.tsv.gz","atac_A5B5"="/home/anjali5/scratch/graham_atac_analysis/atac_A5B5/fragments.tsv.gz","atac_A6B6"="/home/anjali5/scratch/graham_atac_analysis/atac_A6B6/fragments.tsv.gz","atac_A7B7"="/home/anjali5/scratch/graham_atac_analysis/atac_A7B7/fragments.tsv.gz","atac_A8B8"="/home/anjali5/scratch/graham_atac_analysis/atac_A8B8/fragments.tsv.gz","atac_A9B9"="/home/anjali5/scratch/graham_atac_analysis/atac_A9B9/fragments.tsv.gz","atac_A10B10"="/home/anjali5/scratch/graham_atac_analysis/atac_A10B10/fragments.tsv.gz","atac_A11B11"="/home/anjali5/scratch/graham_atac_analysis/atac_A11B11/fragments.tsv.gz","atac_A12B12"="/home/anjali5/scratch/graham_atac_analysis/atac_A12B12/fragments.tsv.gz","atac_A13B13"="/home/anjali5/scratch/graham_atac_analysis/atac_A13B13/fragments.tsv.gz","atac_A14B14"="/home/anjali5/scratch/graham_atac_analysis/atac_A14B14/fragments.tsv.gz","atac_A15B15"="/home/anjali5/scratch/graham_atac_analysis/atac_A15B15/fragments.tsv.gz","atac_A16B16"="/home/anjali5/scratch/graham_atac_analysis/atac_A16B16/fragments.tsv.gz","atac_A17B17"="/home/anjali5/scratch/graham_atac_analysis/atac_A17B17/fragments.tsv.gz","atac_A18B18"="/home/anjali5/scratch/graham_atac_analysis/atac_A18B18/fragments.tsv.gz","atac_A19B19"="/home/anjali5/scratch/graham_atac_analysis/atac_A19B19/fragments.tsv.gz","atac_A20B20"="/home/anjali5/scratch/graham_atac_analysis/atac_A20B20/fragments.tsv.gz","atac_A21B21"="/home/anjali5/scratch/graham_atac_analysis/atac_A21B21/fragments.tsv.gz","atac_A22B22"="/home/anjali5/scratch/graham_atac_analysis/atac_A22B22/fragments.tsv.gz","atac_A23B23"="/home/anjali5/scratch/graham_atac_analysis/atac_A23B23/fragments.tsv.gz","atac_A24B24"="/home/anjali5/scratch/graham_atac_analysis/atac_A24B24/fragments.tsv.gz","atac_A25B25"="/home/anjali5/scratch/graham_atac_analysis/atac_A25B25/fragments.tsv.gz","atac_A26B26"="/home/anjali5/scratch/graham_atac_analysis/atac_A26B26/fragments.tsv.gz","atac_A27B27"="/home/anjali5/scratch/graham_atac_analysis/atac_A27B27/fragments.tsv.gz","atac_A28B28"="/home/anjali5/scratch/graham_atac_analysis/atac_A28B28/fragments.tsv.gz","atac_A29B29"="/home/anjali5/scratch/graham_atac_analysis/atac_A29B29/fragments.tsv.gz","atac_A30B30"="/home/anjali5/scratch/graham_atac_analysis/atac_A30B30/fragments.tsv.gz","atac_A31B31"="/home/anjali5/scratch/graham_atac_analysis/atac_A31B31/fragments.tsv.gz","atac_A32B32"="/home/anjali5/scratch/graham_atac_analysis/atac_A32B32/fragments.tsv.gz","atac_A33B33"="/home/anjali5/scratch/graham_atac_analysis/atac_A33B33/fragments.tsv.gz","atac_A34B34"="/home/anjali5/scratch/graham_atac_analysis/atac_A34B34/fragments.tsv.gz","atac_A35B35"="/home/anjali5/scratch/graham_atac_analysis/atac_A35B35/fragments.tsv.gz","atac_A36B36"="/home/anjali5/scratch/graham_atac_analysis/atac_A36B36/fragments.tsv.gz","atac_A37B37"="/home/anjali5/scratch/graham_atac_analysis/atac_A37B37/fragments.tsv.gz","atac_A38B38"="/home/anjali5/scratch/graham_atac_analysis/atac_A38B38/fragments.tsv.gz","atac_A39B39"="/home/anjali5/scratch/graham_atac_analysis/atac_A39B39/fragments.tsv.gz","atac_A40B40"="/home/anjali5/scratch/graham_atac_analysis/atac_A40B40/fragments.tsv.gz","atac_A41B41"="/home/anjali5/scratch/graham_atac_analysis/atac_A41B41/fragments.tsv.gz","atac_A42B42"="/home/anjali5/scratch/graham_atac_analysis/63mixed10_MPS12346151/outs/fragments.tsv.gz")
ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = names(inputFiles),minTSS = 3,minFrags = 1000,addTileMat = TRUE,addGeneScoreMat = TRUE, excludeChr = "chrM", subThreading = FALSE, threads = 1)
ArrowFiles
#doubScores <- addDoubletScores(input = ArrowFiles, k = 10,knnMethod = "UMAP",LSIMethod = 1)
proj2 <- ArchRProject(ArrowFiles = ArrowFiles,  outputDirectory = "ArchR_proj0",copyArrows = FALSE)
#proj <- loadArchRProject(path = "~/scratch/ArchR_outputs/ArchR_proj1", force = FALSE,showLogo = FALSE)
df <- data.frame(proj@cellColData)
tmp <- do.call('rbind', strsplit(as.character(rownames(df)),'#',fixed=TRUE))
tmp <- as.data.frame(tmp)
df$barcode <- tmp$X2

get_barcodes <- function(val){
  tmp1 <- paste("atac",val, sep="_")
  tmp2 <- dplyr::filter(df, Sample==tmp1)
  tmp1 <- data.frame(tmp2$barcode)
  tmp1$is_cell <- "1"
  names(tmp1) <- c("barcode", "is_cell")
  return(tmp1)
}

list=c("A1B1", "A2B2", "A3B3", "A4B4", "A5B5", "A6B6", "A7B7", "A8B8", "A9B9", "A10B10", "A11B11", "A12B12", "A13B13", "A14B14", "A15B15", "A16B16", "A17B17", "A18B18", "A19B19", "A20B20", "A21B21", "A22B22", "A23B23", "A24B24", "A25B25", "A26B26", "A27B27", "A28B28", "A29B29", "A30B30", "A31B31", "A32B32", "A33B33", "A34B34", "A35B35", "A36B36", "A37B37", "A38B38", "A39B39", "A40B40", "A41B41", "A42B42")

for(i in list){
  assign(paste("atac",i,sep="_"), get_barcodes(i))
  print(head(paste("atac",i,sep="_")))
  write.csv(get(paste("atac",i,sep="_")),paste0("~/scratch/ArchR_outputs/barcodes","atac",i), row.names = FALSE, quote=FALSE)
}






