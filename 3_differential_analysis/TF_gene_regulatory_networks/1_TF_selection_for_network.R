
library(Seurat)
library(tidyverse)
library("ggsci")
library(RColorBrewer)


#load snRNA-seq imputed gene expression matrices for ExN and Mic cell-types
load("~/scratch/TF_gene_exp_selection/motif_geneExp.rds") 

exn <- NormalizeData(exn)
Idents(exn) <- "SubClusters"
ex <- subset(exn,idents="ExN1")

#select top most significant matches from families with multiple motifs and filter by those present in the motif matrix and corresponding gene-exp matrix
#bgzip,zf,BHLH,HTH
breaksList = seq(-1, 1, by = 0.5)

Idents(ex) <- "orig.ident"
features1 <- c("BATF","ATF3","BACH1","FOSL2","BACH2")
features2 <- c("CTCFL","WT1","BCL6","ZNF263","EGR1","CTCF") #,"WT1")
features3 <- c("ATOH1","OLIG2","ASCL1","NEUROD1") #"NEUROG2", "OLIG2")
features4 <- c("RFX6","RFX5","RFX2","RFX1","RFX3")

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
  #print(tmp %>% group_by(features) %>% summarise(mean=median(feature)))
  tmp$features <- factor(tmp$features, levels = unique(tmp$features))
  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/", motifs, "_.pdf"), width=6, height=2)
  print(ggplot(tmp, aes(x=feature,fill=features)) + geom_density(alpha=0.3) + scale_x_log10() + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_fill_brewer(palette = "YlOrRd"))
  dev.off()
}
#########################microglia##################
mic <- NormalizeData(mic)
Idents(mic) <- "SubClusters"
mic2 <- subset(mic,idents="Mic2")

features1 <- c("KLF14","KLF5","KLF9","SP1","KLF3") 
features2 <- c("IRF4","IRF1","IRF2") #"IRF8","ISRE", 
features3 <- c("SPIB","ETV1","ETS1","SPI1","FLI1")
features4 <- c("ATF1","ATF2","CREM","ATF7","JUND")

total_df <- data.frame()
for(features in c("features1", "features2", "features3", "features4")){
  print(features)
  feature_df <- data.frame()
  for(feature in get(features)){
    print(feature)
    #pdf(paste0("/home/anjali5/scratch/TF_gene_exp_selection/","ExN1",features,"Filtered_Vlnplot_TF_selection_for_networks.pdf"))
    VlnPlot(mic2, features = feature) -> d
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
  #print(tmp %>% group_by(features) %>% summarise(mean=median(feature)))
  tmp$features <- factor(tmp$features, levels = unique(tmp$features))
  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/","Mic2", motifs, "_.pdf"), width=6, height=2)
  print(ggplot(tmp, aes(x=feature,fill=features)) + geom_density(alpha=0.3) + scale_x_log10() + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_fill_brewer(palette = "YlOrRd"))
  dev.off()
}
#################################
#ast <- NormalizeData(ast)
#Idents(ast) <- "SubClusters"
#ast3 <- subset(ast,idents="Ast3")
#features1 <- c("AR","NR3C1","NR3C2") #,"PGR") #,"ARE","PR")#
#features2 <- c("CTCFL","CTCF","BCL6")
#features3 <- c("BATF","ATF3","FOSL2","JUN", "BACH2")

#total_df <- data.frame()
#for(features in c("features1", "features2", "features3")){
#  print(features)
#  feature_df <- data.frame()
#  for(feature in get(features)){
#    print(feature)
#    #pdf(paste0("/home/anjali5/scratch/TF_gene_exp_selection/","ExN1",features,"Filtered_Vlnplot_TF_selection_for_networks.pdf"))
#    VlnPlot(ast3, features = feature) -> d
#    colnames(d$data) <- c("feature","ident")
#    d$data[["features"]] <- feature
#    feature_df <- rbind(feature_df, d$data)
#    print(dim(feature_df))
#  }
#  feature_df[["family"]] <- features
#  total_df <- rbind(total_df, feature_df)
#  print(dim(total_df))
#}
#for(motifs in unique(total_df$family)){
#  tmp <- total_df[total_df$family %in% motifs,]
#  #print(tmp %>% group_by(features) %>% summarise(mean=median(feature)))
#  tmp$features <- factor(tmp$features, levels = unique(tmp$features))
#  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/","Ast3", motifs, "_.pdf"), width=6, height=2)
#  print(ggplot(tmp, aes(x=feature,fill=features)) + geom_density(alpha=0.3) + scale_x_log10() + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_fill_brewer(palette = "YlOrRd"))
#  dev.off()
#}

#########################################
#https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


########split violin plots for chromVar TF z-score (Case-Control) for selection of network TFs##############
library(ArchR)
proj <- loadArchRProject("~/scratch/ArchR_backup/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")
motifPositions <- getPositions(proj, name="Motif")
proj2 <- proj[proj$SubClusters %in% "ExN1",]

features1 <- c("FOSL2","BACH2")
features2 <- c("EGR1","CTCF") 
features3 <- c("ASCL1","NEUROD1")
features4 <- c("RFX1","RFX3")

for(features in c("features1","features2","features3","features4")){
  feature_df <- data.frame()
  for(feature in get(features)){
    print(feature)
    markerMotif <- unlist(lapply(feature, function(x) grep(x, names(motifPositions), value = TRUE)))[1]
    p <- plotGroups(proj2, name=paste0("z:",markerMotif), groupBy="condition", colorBy="MotifMatrix")
    #print(p)
    p$data[["feature"]] <- feature
    feature_df <- rbind(feature_df, p$data)
    #print(feature_df)
  }
  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/","ExN1", features,"nop_chromVar.pdf"), width=5, height=4)
  print(ggplot(feature_df, aes(feature, y, fill = x)) + geom_split_violin() + scale_fill_manual(values=c("tan","salmon4")) + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()))
  dev.off()
}
# ggpubr::stat_compare_means(method="wilcox.test",label =NULL) +
#############
features1 <- c("SP1","KLF3") 
features2 <- c("IRF1","IRF2") 
features3 <- c("SPI1","FLI1")
features4 <- c("ATF7","JUND")


proj2 <- proj[proj$SubClusters %in% "Mic2",]

for(features in c("features1","features2","features3","features4")){
  feature_df <- data.frame()
  for(feature in get(features)){
    print(feature)
    markerMotif <- unlist(lapply(feature, function(x) grep(x, names(motifPositions), value = TRUE)))[1]
    p <- plotGroups(proj2, name=paste0("z:",markerMotif), groupBy="condition", colorBy="MotifMatrix")
    #print(p)
    p$data[["feature"]] <- feature
    feature_df <- rbind(feature_df, p$data)
    #print(feature_df)
  }
  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/","Mic2",features,"nop_chromVar.pdf"), width=5, height=4)
  print(ggplot(feature_df, aes(feature, y, fill = x)) + geom_split_violin() + scale_fill_manual(values=c("tan","salmon4")) + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()))
  dev.off()
}

################
#features1 <- c("AR","NR3C1","NR3C2") 
#features2 <- c("CTCF","BCL6")
#features3 <- c("JUN", "BACH2")

#proj2 <- proj[proj$SubClusters %in% "Ast3",]

#for(features in c("features1","features2","features3","features4")){
#  feature_df <- data.frame()
#  for(feature in get(features)){
#    print(feature)
#    markerMotif <- unlist(lapply(feature, function(x) grep(x, names(motifPositions), value = TRUE)))[1]
#    p <- plotGroups(proj2, name=paste0("z:",markerMotif), groupBy="condition", colorBy="MotifMatrix")
#    #print(p)
#    p$data[["feature"]] <- feature
#    feature_df <- rbind(feature_df, p$data)
    #print(feature_df)
#  }
#  pdf(paste0("~/scratch/TF_gene_exp_selection/Updated_TF_selection/","Ast3", features,"nop_chromVar.pdf"), width=5, height=4)
#  print(ggplot(feature_df, aes(feature, y, fill = x)) + geom_split_violin() + scale_fill_manual(values=c("tan","salmon4")) + theme_classic() + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()))
#  dev.off()
#}

#ggpubr::stat_compare_means(method="wilcox.test", label ="p.format") + 

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

