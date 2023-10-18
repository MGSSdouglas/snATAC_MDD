library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(tidyverse)
library(data.table)
library(patchwork)
library(TFBSTools)
library(JASPAR2020)
library(SeuratObject)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)

#Load result file
atac <- readRDS("~/projects/def-gturecki/anjali5/Finalized_DAR/subcluster/res_freq_subcluster.rds")
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct
#get significant DARs
atac_sub <- atac[atac$p_adj.loc_treatment < 0.05 & abs(atac$logFC_treatment) > log2(1.1) & atac$Greater3pct==TRUE,]
length(unique(atac_sub$gene))

#convert per-cluster DARs split by logFC direction into bed files ready for homer
clusters= unique(atac_sub$cluster_id)
directions = c("up", "down")
for(cluster in clusters){
  print(cluster)
  table_cur <- atac_sub[atac_sub$cluster_id %in% cluster,]
  if(nrow(table_cur) > 0){
    for(direction in directions){
      print(direction)
      if(direction=="up"){
        cur <- dplyr::filter(table_cur, logFC_treatment > 0)
      } else {
        cur <- dplyr::filter(table_cur, logFC_treatment < 0)
      }
      print(dim(cur))
      top.da.peaks <- cur$gene
      if(length(top.da.peaks) >  0){
        #convert string to bed file
        output <- "/home/anjali5/scratch/"
        y <- do.call(rbind.data.frame,strsplit(unique(top.da.peaks), split="-"))
        print(dim(y))
        y$peakID <- paste0("Peak",seq.int(nrow(y)))
        y$unuse <- rep(".",nrow(y))
        y$strand <- rep(0,nrow(y))
        print(dim(y))
        write.table(y, paste0(output,cluster,direction,"_peak.bed"),row.names=FALSE, col.names = FALSE,quote = FALSE,sep="\t")
      }
    }
  }
  
