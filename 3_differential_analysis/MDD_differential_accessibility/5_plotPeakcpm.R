#result <- readRDS("~/Desktop/subcluster_result_limma_allPeaks_corrected_filtered.rds")
#plot E-values of limma-voom from differential analysis for each group

library(readxl)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(tidyverse)

atac2 <- readRDS("~/Desktop/res_freq_broad.rds")
tmp   <- atac2[, grep("Case.frq|Control.frq", colnames(atac2))]
tmp$Greater3pct <- ifelse(tmp$Case.frq  > 0.01 | tmp$Control.frq > 0.01, "TRUE", "FALSE")
atac2$Greater3pct <- tmp$Greater3pct


atac <- readRDS("~/Desktop/res_freq_subcluster.rds")
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct

meta <- read_excel("~/Desktop/psychencode_data/metadata_subjects.xlsx")
meta <- meta[meta$Subject != "atac_A32B32female",]

#conclusions <- unique(atac$conclusion)
#clusters=c("ExN1","ExN2")
#genes=c("chr22-41091360-41091860", "chr2-234446364-234446864")
#####plot list of peaks ########
for (gene in genes){
  print(gene)
  for(cluster in clusters){
    print(cluster)
    #plots_list <- list()
    atac_tmp <- atac[atac$gene==gene & atac$cluster_id == cluster,]
    if(nrow(atac_tmp) > 0){
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
      tmp$interaction <- NULL
      colnames(tmp) <- c("peak","condition")
      pdf(file=paste0("~/Desktop/atac_2023_paper/peaks/",cluster, gene, ".pdf"), onefile = TRUE, width=5, height=7)
      print(ggplot(tmp, aes(condition,peak,fill=condition)) + geom_boxplot() + theme_classic() + scale_fill_manual(values=c('salmon4','tan'))+ theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 90, size = 20),axis.text.y = element_text(angle = 0, size = 20)) + ylab(paste0(gene," ", "(cpm)")))
      dev.off()
    }
  }
}











