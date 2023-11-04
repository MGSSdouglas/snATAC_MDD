#result <- readRDS("~/Desktop/subcluster_result_limma_allPeaks_corrected_filtered.rds")
# This script will plot limma voom E values

library(readxl)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(tidyverse)

atac <- readRDS("~/Desktop/res_freq_subcluster.rds")
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct

meta <- read_excel("~/Desktop/metadata_subjects.xlsx")
meta <- meta[meta$Subject != "atac_A32B32female",]

plots_list <- list()
clusters <- c("ExN1", "Mic2") #, "Ast3")
directions <- c("up","down")
for (cluster in clusters){
  print(cluster)
  for(direction in directions){
    if(direction == "up"){
    atac_tmp <- atac[atac$logFC_treatment > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE & atac$cluster_id==cluster,]
    } else{
    atac_tmp <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE & atac$cluster_id=="Ast3",]
    }
    if(nrow(atac_tmp) > 0){
      print(direction)
      tmp <- atac_tmp[, grep(".cpm|gene", colnames(atac_tmp))]
      rownames(tmp) <- tmp$gene
      tmp$gene <- NULL
      tmp <- t(tmp)
      tmp <- as.data.frame(tmp)
      rownames(tmp) <- sub(".cpm","",rownames(tmp))
      meta <- meta[order(match(meta$Subject,rownames(tmp))),]
      tmp$condition <- meta$Causeofdeath
      tmp$condition <- ifelse(tmp$condition=="Suicide" | tmp$condition=="Case", "case", "control")
      tmp %>% na.omit() ->  tmp
      tmp %>% arrange(condition) -> tmp
      group <- data.frame(tmp$condition)
      rownames(group) <- rownames(tmp)
      tmp$condition <- NULL
      colnames(group) <- "groups"
      #group_sub <- subset(group,groups %in% c("male_case", "male_control","female_case","female_control"))
      group_sub <- subset(group,groups %in% c("case", "control"))
      tmp_sub <- subset(tmp, rownames(tmp) %in% rownames(group_sub))
          chrom <- do.call(rbind,stringr::str_split(colnames(tmp_sub),"-"))
          rownames(chrom) <- colnames(tmp_sub)
          colnames(chrom) <- c("Chromosome", "start", "end")
          chrom <- data.frame(chrom)
          chrom$start <- NULL
          chrom$end <- NULL
          pheatmap::pheatmap(tmp_sub, cluster_cols = F, cluster_rows = F, annotation_row  = group, annotation_col = chrom, show_colnames  = FALSE)
          pdf(file=paste0("~/Desktop/plots_atac/combined_heatmaps_",cluster,direction,".pdf"))
          pheatmap::pheatmap(tmp_sub, scale = "column", cluster_cols = F, cluster_rows = F, annotation_row  = group_sub, gaps_row = as.numeric(cumsum(table(group_sub))),main=paste0(cluster,"-",direction), annotation_col = chrom, show_colnames  = FALSE) #, color=colorRampPalette(c("white", "beige", "brown"),bias=1)(200))
          dev.off()
    }
  }
}



