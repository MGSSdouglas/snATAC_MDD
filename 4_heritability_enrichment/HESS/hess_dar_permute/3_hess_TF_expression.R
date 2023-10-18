
library(Seurat)
library(tidyverse)
library("ggsci")
library(RColorBrewer)

#density plots for hess tfs enriched in ExN1 
features1 <- rev(c("TCF4","RBPJ","TWIST2","NEUROD1","OLIG2","NEUROG2","ATOH1","BHLHA15"))
features2 <- rev(c("BACH2","FOSL2","JUN","JUNB","BACH1","FOS","ATF3" ,"MAFA","BATF"))
features3 <-rev(c("RFX3","RFX1","RFX2","RFX4","RFX5","RFX6"))
features4 <- rev(c("CTCF","EGR1","ZNF341","EGR2","HIC1","WT1","CTCFL","PRDM9","KLF14"))

Idents(ex) <- "orig.ident"
total_df <- data.frame()
total_df <- data.frame()
for(features in c("features1", "features2", "features3", "features4")){
  print(features)
  feature_df <- data.frame()
  for(feature in get(features)){
    print(feature)
    #pdf(paste0("/home/anjali5/scratch/TF_gene_exp_selection/","ExN1",features,"Filtered_Vlnplot_TF_selection_for_networks.pdf"))
    VlnPlot(ex, features = feature) -> d
    colnames(d$data) <- c("feature","ident")
    d$data[["features"]] <- feature
    feature_df <- rbind(feature_df, d$data)
    print(dim(feature_df))
  }
  feature_df[["family"]] <- features
  total_df <- rbind(total_df, feature_df)
  print(dim(total_df))
}
for(motifs in unique(total_df$family)){
  tmp <- total_df[total_df$family %in% motifs,]
  #print(tmp %>% group_by(features) %>% summarise(median=median(feature)) %>% arrange(desc(median)))
  tmp$features <- factor(tmp$features, levels = unique(tmp$features))
  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/HESS_TF_genes/",motifs, ".pdf"), width=6, height=3)
  print(ggplot(tmp, aes(x=feature,fill=features)) + geom_density(alpha=0.3) + scale_x_log10() + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_fill_brewer(palette = "YlOrRd"))
  dev.off()
}

#plot..
Idents(ex) <- "orig.ident"
DotPlot(ex, features = features) + RotatedAxis() -> d
d$data %>% select(features.plot,avg.exp) %>% pivot_wider(names_from=features.plot,values_from=avg.exp) ->d
pdf("~/scratch/bhlh_unscaled.pdf")
pheatmap::pheatmap(d, cluster_rows=F,cluster_cols=F, scale="none",show_rownames=F,cellheight=10,color=colorRampPalette(c("white","peachpuff1", "red4"))(50),fontsize=20)
dev.off()

