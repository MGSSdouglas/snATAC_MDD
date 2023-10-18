
library(Seurat)
library(tidyverse)
library("ggsci")
library(RColorBrewer)


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


