
library(Seurat)
library(tidyverse)
library(reshape2)
library(ggpubr)

##########load MDD snRNA-seq data########################
#load("/home/anjali5/projects/def-gturecki/anjali5/module_scores/rna_module_added.rds")
#load("/home/anjali5/projects/def-gturecki/anjali5/module_scores/rna_module_updated.rds")
snRNA <- readRDS("/home/anjali5/projects/def-cnagy/1_enhanced_harmonized_object.Rds")
Idents(snRNA) <- "Cluster"
exn_deep <- subset(snRNA, idents=c("ExN16_L56"))
mic_rna <- subset(snRNA, idents=c("Mic1"))
###########load TF-target genes for Module score##########
ExN1_genes <- read.csv("/home/anjali5/scratch/TF_gene_exp_selection/Updated_TF_selection/Updated_ExN1_edge_df.csv")
#Ast3_genes <- read.csv("/home/anjali5/scratch/TF_gene_exp_selection/Updated_TF_selection/Updated_Ast3_edge_df.csv")
Mic2_genes <- read.csv("/home/anjali5/scratch/TF_gene_exp_selection/Updated_TF_selection/Updated_Mic2_edge_df.csv")

exn_deep <- AddModuleScore(exn_deep, features=list("Tf1"=ExN1_genes$to[ExN1_genes$from %in% "Bach2.bZIP"],
                                                     "Tf2"=ExN1_genes$to[ExN1_genes$from %in% "CTCF.Zf"],
                                                     "Tf3"=ExN1_genes$to[ExN1_genes$from %in% "Ascl1.bHLH"],
                                                     "Tf4"=ExN1_genes$to[ExN1_genes$from %in% "X.box.HTH"],
                                                      "total"=ExN1_genes$to),
                                                     name = "Tfactor")
mic_rna <- AddModuleScore(mic_rna, features=list("Tf1"=Mic2_genes$to[Mic2_genes$from %in% "PU.1.ETS"],
                                                    "Tf2"=Mic2_genes$to[Mic2_genes$from %in% "IRF1.IRF"],
                                                   "Tf3"=Mic2_genes$to[Mic2_genes$from %in% "Sp1.Zf"],
                                                   "Tf4"=Mic2_genes$to[Mic2_genes$from %in% "Atf7.bZIP"],
                                                 "total"=Mic2_genes$to),
                                                     name = "Tfactor")

#ast_rna <- AddModuleScore(ast_rna, features=list("Tf1"=Ast3_genes$to[Ast3_genes$from %in% "GRE.NR..IR3"],
#                                                 "Tf2"=Ast3_genes$to[Ast3_genes$from %in% "AR.halfsite.NR"],
#                                                 "Tf3"=Ast3_genes$to[Ast3_genes$from %in% "Bach2.bZIP"],
#                                                 "Tf4"=Ast3_genes$to[Ast3_genes$from %in% "CTCF.Zf"],
#                                                 "total"=Ast3_genes$to, "Tf5"=Ast3_genes$to[Ast3_genes$from %in% c("GRE.NR..IR3","AR.halfsite.NR")],
#                                                 "Tf6"=Ast3_genes$to[Ast3_genes$from %in% c("CTCF.Zf","Bach2.bZIP")]),name = "Tfactor")
###############################Geom split violin###################################################
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

################Prepare to plot#################
Idents(exn_deep) <- "Condition"
tmpdf <- list()
for(tf in c("Tfactor1","Tfactor2","Tfactor3","Tfactor4","Tfactor5")){
  p <- VlnPlot(exn_deep,features = tf)
  tmpdf[[tf]] <- p$data
}
total_df <- data.frame(tmpdf)

total_df %>% select(Tfactor1.ident,Tfactor1.Tfactor1,Tfactor2.Tfactor2,Tfactor3.Tfactor3,Tfactor4.Tfactor4,Tfactor5.Tfactor5) %>% reshape2::melt() -> plot_df
plot_df$variable <- ifelse(plot_df$variable=="Tfactor1.Tfactor1", "BACH2", ifelse(plot_df$variable=="Tfactor2.Tfactor2", "CTCF", ifelse(plot_df$variable=="Tfactor3.Tfactor3", "ASCL1", ifelse(plot_df$variable=="Tfactor4.Tfactor4","RFX3","ALL"))))
#plot_df %>% mutate(variable = fct_relevel(variable,"BACH2","CTCF","ASCL1","RFX3","ALL")) %>%ggplot(aes(x=variable,y=value,fill=Tfactor1.ident)) + geom_split_violin() + stat_compare_means(method="wilcox.test") +theme_classic() + scale_fill_manual(values=c("salmon4","tan")) + theme(text = element_text(size = 20),axis.text = element_text(size = 20))

for(var in unique(plot_df$variable)){
  print(var)
  tmp <- plot_df[plot_df$variable %in% var,]
  pdf(paste0("/home/anjali5/scratch/TF_gene_exp_selection/Updated_TF_selection/", "ExN_deep_", var, ".pdf"),width=4, height=4)
  print(ggplot(tmp, aes(x=variable,y=value,fill=Tfactor1.ident)) + geom_split_violin() + stat_compare_means(method="wilcox.test", label="p.format") +theme_classic() + scale_fill_manual(values=c("salmon4","tan")) + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_y_continuous(breaks=seq(-0.5,0.5,by=0.3)))
  dev.off()
}

################repeat for microglia####################
Idents(mic_rna) <- "Condition"
tmpdf <- list()
for(tf in c("Tfactor1","Tfactor2","Tfactor3","Tfactor4","Tfactor5")){
  p <- VlnPlot(mic_rna,features = tf)
  tmpdf[[tf]] <- p$data
}
total_df <- data.frame(tmpdf)

total_df %>% select(Tfactor1.ident,Tfactor1.Tfactor1,Tfactor2.Tfactor2,Tfactor3.Tfactor3,Tfactor4.Tfactor4,Tfactor5.Tfactor5) %>% reshape2::melt() -> plot_df
plot_df$variable <- ifelse(plot_df$variable=="Tfactor1.Tfactor1", "SPI1", ifelse(plot_df$variable=="Tfactor2.Tfactor2", "IRF1", ifelse(plot_df$variable=="Tfactor3.Tfactor3", "SP1", ifelse(plot_df$variable=="Tfactor4.Tfactor4","ATF7","ALL"))))

for(var in unique(plot_df$variable)){
  print(var)
  tmp <- plot_df[plot_df$variable %in% var,]
  pdf(paste0("/home/anjali5/scratch/TF_gene_exp_selection/Updated_TF_selection/","Mic2_", var, ".pdf"),width=4, height=4)
  print(ggplot(tmp, aes(x=variable,y=value,fill=Tfactor1.ident)) + geom_split_violin() + stat_compare_means(method="wilcox.test", label="p.format") +theme_classic() + scale_fill_manual(values=c("salmon4","tan")) + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_y_continuous(breaks=seq(-0.5,0.5,by=0.3)))
  dev.off()
}
################
#Idents(ast_rna) <- "Condition"
##tmpdf <- list()
#for(tf in c("Tfactor1","Tfactor2","Tfactor3","Tfactor4","Tfactor5","Tfactor6","Tfactor7")){
#  p <- VlnPlot(ast_rna,features = tf,split.by = "Cluster")
#  tmpdf[[tf]] <- p$data
#}
#total_df <- data.frame(tmpdf)

#total_df %>% select(Tfactor1.split,Tfactor1.ident,Tfactor1.Tfactor1,Tfactor2.Tfactor2,Tfactor3.Tfactor3,Tfactor4.Tfactor4,Tfactor5.Tfactor5,Tfactor6.Tfactor6,Tfactor7.Tfactor7) %>% reshape2::melt() -> plot_df
#plot_df$variable <- ifelse(plot_df$variable=="Tfactor1.Tfactor1", "NR3C1", ifelse(plot_df$variable=="Tfactor2.Tfactor2", "AR", ifelse(plot_df$variable=="Tfactor3.Tfactor3", "BACH2", ifelse(plot_df$variable=="Tfactor4.Tfactor4","CTCF", ifelse(plot_df$variable=="Tfactor5.Tfactor5","ALL", ifelse(plot_df$variable=="Tfactor6.Tfactor6","Up", "Down"))))))

#for(var in unique(plot_df$variable)){
#  print(var)
#  tmp <- plot_df[plot_df$variable %in% var,]
#  print(tmp %>% group_by(Tfactor1.ident) %>% summarise(mean=mean(value)))
#  print(tmp %>% group_by(Tfactor1.ident) %>% summarise(mean=median(value)))
#  pdf(paste0("/home/anjali5/scratch/TF_gene_exp_selection/Updated_TF_selection/","Ast_", var, ".pdf"),width=4, height=4)
#  print(ggplot(tmp, aes(x=variable,y=value,fill=Tfactor1.ident)) + geom_split_violin() + stat_compare_means(method="wilcox.test", label="p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.1, Inf), symbols = c( "***", "**", "*","ns"))) +theme_classic() + scale_fill_manual(values=c("salmon4","tan")) + theme(text = element_text(size = 20),axis.text = element_text(size = 20),legend.title=element_blank()) + scale_y_continuous(breaks=seq(-0.5,0.8,by=0.3))) #+ facet_wrap(~Tfactor1.split))
#  dev.off()
#}

#symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))


