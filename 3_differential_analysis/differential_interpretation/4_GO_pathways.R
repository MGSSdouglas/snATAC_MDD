library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GeneOverlap)

#load atac data
atac <- readRDS("~/Desktop/res_freq_subcluster.rds")
#atac$cluster <- plyr::mapvalues(atac$cluster_id, from = unique(atac$cluster_id), to = c("Ast", "Ast", "Ast", "Ast", "Mic", "ExN", "ExN","ExN", "ExN","ExN", "ExN","ExN", "ExN", "ExN", "ExN", "InN", "InN", "InN", "InN", "Mic", "Mic", "Oli", "Oli", "Oli", "Oli", "Oli", "OPC", "OPC", "OPC"))
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE,] #| abs(atac$logFC_male.trt) > log2(1.1) & atac$p_adj.loc_male.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_female.trt) > log2(1.1) & atac$p_adj.loc_female.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_interaction) > log2(1.1) & atac$p_adj.loc_interaction < 0.2 & atac$Greater3pct==TRUE,]

#get p2g linkages
sub_p2g <- readRDS("~/Desktop/Subtype_p2g_redone.rds")
p2g <- sub_p2g %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45 & FDR < 1e-4)

##load background genes
gc <- readRDS("~/Desktop/IterativelyRemovedPeaks_sub.rds")
#gc$peak <- paste(seqnames(gc), ranges(gc), sep="-")
#loop over every cluster and get pathways in atac
for(cluster in unique(atac_sub$cluster_id)){
  print(cluster)
  atac_tmp <- atac_sub$gene[atac_sub$cluster_id %in% cluster]
  atac_genes <- p2g$geneName[p2g$peakName %in% atac_tmp]
  print(length(atac_genes))
  res_atac <- enrichGO(gene = atac_genes, universe= unique(gc$nearestGene), OrgDb = org.Hs.eg.db, ont = "ALL",  keyType = 'SYMBOL', pAdjustMethod = "BH",pvalueCutoff = 0.05,minGSSize = 5)
  res_atac
  assign(paste0("atac",cluster), res_atac)
}

#> saveRDS(atacExN1, "~/Desktop/atac_2023_paper/ExN1_pathways.rds")
#> saveRDS(atacMic2, "~/Desktop/atac_2023_paper/Mic2_pathways.rds")

#plot ExN1 and Mic2 pathway networks
library(enrichplot)
edox <- setReadable(atacMic2, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
#p1 <- treeplot(edox2, fontsize = 8)

pdf("~/Desktop/atac_2023_paper/mic2_pathway_emapplot.pdf", width = 8)
emapplot(edox2,cex_label_category = 0.7,min_edge = 0.1,legend_n = 3) + viridis::scale_fill_viridis(option = "A")
dev.off()

pdf("~/Desktop/atac_2023_paper/exn1_pathway_emapplot.pdf", width = 8)
emapplot(edox2,cex_label_category = 0.7,min_edge = 0.1,legend_n = 3) + viridis::scale_fill_viridis(option = "A")
dev.off()


#get pathways split by direction
#loop over every cluster and get pathways in atac 
for(cluster in unique(atac_sub$cluster_id)){
  print(cluster)
  atac_1 <- atac_sub[atac_sub$cluster_id %in% cluster,]
  #atac
  for(direction in c("up","down")){
    if(direction == "up"){
      atac_peaks <- atac_1$gene[atac_1$logFC_treatment >0]
    }else{
      atac_peaks <- atac_1$gene[atac_1$logFC_treatment <0]
    }
    atac_genes <- p2g$geneName[p2g$peakName %in% atac_peaks]
    print(length(atac_genes))
    if(length(atac_genes) >0){
      #eg = bitr(atac_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      res_atac <- enrichGO(gene = atac_genes, universe = gc$nearestGene,OrgDb = org.Hs.eg.db, ont = "BP",  keyType = 'SYMBOL', pAdjustMethod = "BH",pvalueCutoff = 0.2, minGSSize = 5,qvalueCutoff = 1)
      res_atac
      assign(paste0("atac",cluster, direction), res_atac)
    }
  }
}
for(i in unique(c(ls()))){
  print(i)
  print(get(i))
}
total_df <- data.frame()
for(obj in ls()){
  print(obj)
  print(get(obj))
  tmp <- get(obj)@result
  if(nrow(tmp) > 0){
    tmp$cluster <-  as.character(obj)
    total_df=rbind(total_df, tmp)
  }
}

#write to output
write.csv(total_df, "~/scratch/cluster_pathways.csv")

#repeat for broad
atac2 <- readRDS("~/Desktop/res_freq_broad.rds")
tmp   <- atac2[, grep("Case.frq|Control.frq", colnames(atac2))]
tmp$Greater3pct <- ifelse(tmp$Case.frq  > 0.01 | tmp$Control.frq > 0.01, "TRUE", "FALSE")
atac2$Greater3pct <- tmp$Greater3pct
atac_broad <- atac2[abs(atac2$logFC_treatment) > log2(1.1) & atac2$p_adj.loc_treatment < 0.05 & atac2$Greater3pct==TRUE,] #| abs(atac$logFC_male.trt) > log2(1.1) & atac$p_adj.loc_male.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_female.trt) > log2(1.1) & atac$p_adj.loc_female.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_interaction) > log2(1.1) & atac$p_adj.loc_interaction < 0.2 & atac$Greater3pct==TRUE,]

broad_p2g <- readRDS("~/Desktop/broad_p2g.rds")
broad_p2g$peakName <- sub(':', '-', broad_p2g$peakName)
broad_p2g_peaks  <- broad_p2g %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45 & FDR < 1e-4)
broad_p2g_peaks$EmpPval <- NULL
broad_p2g_peaks$EmpFDR <- NULL
p2g <- rbind(broad_p2g_peaks)
##
#loop over every cluster and get pathways in atac and rna
for(cluster in unique(atac_broad$cluster_id)){
  print(cluster)
  atac_1 <- atac_broad[atac_broad$cluster_id %in% cluster,]
  #atac
  for(direction in c("up","down")){
    if(direction == "up"){
      atac_peaks <- atac_1$gene[atac_1$logFC_treatment >0]
    }else{
      atac_peaks <- atac_1$gene[atac_1$logFC_treatment <0]
    }
    atac_genes <- p2g$geneName[p2g$peakName %in% atac_peaks]
    print(length(atac_genes))
    if(length(atac_genes) >0){
      res_atac <- enrichGO(gene = atac_genes, universe= unique(gc$nearestGene), OrgDb = org.Hs.eg.db, ont = "ALL",  keyType = 'SYMBOL', pAdjustMethod = "BH",pvalueCutoff = 0.2, minGSSize = 5,qvalueCutoff = 1)
      res_atac
      assign(paste0("atac",cluster,direction), res_atac)
    }
  }
}

total_df <- data.frame()
for(obj in ls()){
  print(obj)
  print(get(obj))
  tmp <- get(obj)@result
  if(nrow(tmp) > 0){
    tmp$cluster <-  as.character(obj)
    total_df=rbind(total_df, tmp)
  }
}

#write to output
write.csv(total_df, "~/scratch/broad_pathways.csv")


#plot top microglia pathways at broad and cluster levels
###for microglia 
rm(list=setdiff(ls(),c("atacMicdown","atacMicup", "atacMic2down","atacMic1down", "atacMic1up")))

total_df <- data.frame()
for(obj in c("atacMicdown", "atacMicup", "atacMic1down","atacMic1up", "atacMic2down")){
  print(obj)
  name=get(obj)@result$Description[get(obj)@result$p.adjust < 0.2][1:10]
  pvalue=-log10(get(obj)@result$pvalue[get(obj)@result$p.adjust < 0.2][1:10])
  df <- data.frame("Pathway"=name, "Pvalue"=pvalue,"type"=paste(as.character(obj)))
  total_df <- rbind(total_df,df)
}

plot_df <- total_df %>% pivot_wider(names_from = "Pathway", values_from="Pvalue") %>% column_to_rownames("type") 
#plot_df <- plot_df[match(c("atacMicdown", "atacMicup", "atacMic1down","atacMic1up", "atacMic2down"),rownames(plot_df)),]
anno_df <- data.frame("celltype"=c("Mic","Mic" ,"Mic1", "Mic1", "Mic2"),"direction"=c("Down","Up","Down","Up", "Down"))
rownames(anno_df) <- rownames(plot_df)
pal <- ArchR::ArchRPalettes$whitePurple
annoCol<-list(celltype=c("Mic"="cadetblue","Mic1" = "#1792AA", "Mic2"="#83D6E6"),direction=c("Down"="darkslategrey","Up"="goldenrod2"))
pdf("~/scratch/Mic_pathway_heatmap.pdf",width = 25)

library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
pheatmap::pheatmap(plot_df,annotation_row = anno_df,annotation_colors = annoCol, cluster_cols = F, cluster_rows = F, color=colorRampPalette(pal)(256),cellwidth = 20,cellheight = 15,show_rownames = F,fontsize = 15,annotation_names_col=F,annotation_names_row=F)
dev.off()


######marker cCRE pathways######

library(clusterProfiler)
library(org.Hs.eg.db)

sub_p2g_full <- readRDS("~/Desktop/Subtype_p2g_redone.rds")
sub_p2g  <- sub_p2g_full %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45 & FDR < 1e-04)

path="~/Desktop/List_of_SubClusters_Marker_Peaks/"
file_list <- list.files(path=path)

for(file_number in 1:length(file_list)){
  #print(paste0(file_list[file_number]))
  #print(i)
  #if(paste0(file_list[file_number]) %in% "ExN1_markerpeaks.tsv"){
  print(paste0(file_list[file_number]))
  tmp_cluster <- read.csv(paste0(path,file_list[file_number]), sep="\t",header=TRUE)
  tmp_cluster$gene <- paste(tmp_cluster$seqnames, tmp_cluster$start, tmp_cluster$end, sep="-")
  tmp_cluster <- tmp_cluster[order(tmp_cluster$FDR),]
    #print(head(tmp_cluster))
  res <- enrichGO(gene = sub_p2g$geneName[sub_p2g$peakName %in% tmp_cluster$gene], OrgDb = org.Hs.eg.db, ont = "ALL",  keyType = 'SYMBOL', pAdjustMethod = "BH",pvalueCutoff = 0.05, minGSSize = 5)
  }



