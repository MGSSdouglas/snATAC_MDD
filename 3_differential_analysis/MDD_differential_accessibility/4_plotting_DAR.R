library(ggplot2)
library(ggrepel)
library(ggsci)
library(tidyverse)
library(org.Hs.eg.db)
library(GenomicRanges)

#Broad DAR
atac2 <- readRDS("~/Desktop/res_freq_broad.rds")
tmp   <- atac2[, grep("Case.frq|Control.frq", colnames(atac2))]
tmp$Greater3pct <- ifelse(tmp$Case.frq  > 0.01 | tmp$Control.frq > 0.01, "TRUE", "FALSE")
atac2$Greater3pct <- tmp$Greater3pct
atac_broad <- atac2[abs(atac2$logFC_treatment) > log2(1.1) & atac2$p_adj.loc_treatment < 0.05 & atac2$Greater3pct==TRUE,] #| abs(atac$logFC_male.trt) > log2(1.1) & atac$p_adj.loc_male.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_female.trt) > log2(1.1) & atac$p_adj.loc_female.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_interaction) > log2(1.1) & atac$p_adj.loc_interaction < 0.2 & atac$Greater3pct==TRUE,]

##peak metadata##
gc <- readRDS("~/Desktop/IterativelyRemovedPeaks.rds")
gc$peaks <- paste(seqnames(gc), ranges(gc),sep="-")
gc2 <- gc[match(atac_broad$gene, gc$peaks),]
gc2 <- gc2[gc2$peaks %in% atac_broad$gene]
atac_broad$distTSS <- gc2$distToTSS
atac_broad$distToTSS <- ifelse(atac_broad$distTSS < 2000, "<2kbp", ">2kbp")
atac_broad$logFC <- atac_broad$logFC_treatment
pdf("~/Desktop/atac_2023_paper/broad_dar_plot.pdf",width=8)
ggplot(atac_broad, aes(x=cluster_id, y=logFC)) + 
  geom_jitter(aes(colour = distToTSS), size = 1, height = 0.01) +
  geom_segment(aes(x=cluster_id,xend=cluster_id,y=min(logFC),yend=max(logFC)),linetype="dashed", size=0.1) + 
  geom_hline(yintercept = 0, linetype="dotted") + scale_color_manual(values = c("<2kbp" = "#F21A00", ">2kbp" = "#3B9AB2")) + 
  theme_classic() +theme(axis.text = element_text(size = 20,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0), axis.title.x = element_text(size = 20, angle = 0, hjust = .5, vjust = 0), axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),legend.text=element_text(size=20),legend.title=element_text(size=20))
dev.off()

#cluster DAR
atac <- readRDS("~/Desktop/res_freq_subcluster.rds")
atac$cluster_id <- plyr::mapvalues(atac$cluster_id, from = unique(atac$cluster_id), to = c("Ast1", "Ast2", "Ast3", "Ast4", "End1", "ExN1", "ExN11","ExN10", "ExN3","ExN2", "ExN12","ExN8", "ExN4", "ExN9", "ExN5", "InNLAMP5", "InNPV", "InNSST", "InNVIP", "Mic1", "Mic2", "Oli2", "Oli4", "Oli5", "Oli6", "Oli7", "OPC1", "OPC3", "OPC4"))
atac$cluster <- plyr::mapvalues(atac$cluster_id, from = unique(atac$cluster_id), to = c("Ast", "Ast", "Ast", "Ast", "End", "ExN", "ExN","ExN", "ExN","ExN", "ExN","ExN", "ExN", "ExN", "ExN", "InN", "InN", "InN", "InN", "Mic", "Mic", "Oli", "Oli", "Oli", "Oli", "Oli", "OPC", "OPC", "OPC"))
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE,] #| abs(atac$logFC_male.trt) > log2(1.1) & atac$p_adj.loc_male.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_female.trt) > log2(1.1) & atac$p_adj.loc_female.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_interaction) > log2(1.1) & atac$p_adj.loc_interaction < 0.2 & atac$Greater3pct==TRUE,]

##peak metadata##
gc <- readRDS("~/Desktop/IterativelyRemovedPeaks_sub.rds")
gc$peaks <- paste(seqnames(gc), ranges(gc),sep="-")
gc2 <- gc[match(atac_sub$gene, gc$peaks),]
gc2 <- gc2[gc2$peaks %in% atac_sub$gene]
atac_sub$distTSS <- gc2$distToTSS
atac_sub$distToTSS <- ifelse(atac_sub$distTSS < 2000, "<2kbp", ">2kbp")
atac_sub$logFC <- atac_sub$logFC_treatment
atac_sub %>% arrange(factor(cluster, levels = c("Ast","End","Mic","ExN","InN","Oli","OPC"))) -> atac_sub
atac_sub$cluster <- factor(atac_sub$cluster, levels = unique(atac_sub$cluster))
pdf("~/Desktop/atac_2023_paper/Fig2/sub_dar_plot.pdf",width=18,height=8)
#pdf("~/Desktop/tmp.pdf",width=18,height=8)
ggplot(atac_sub, aes(x=cluster_id, y=logFC)) + 
  geom_jitter(aes(colour = distToTSS), size = 1, height = 0.01) +
  geom_segment(aes(x=cluster_id,xend=cluster_id,y=min(logFC),yend=max(logFC)),linetype="dashed", size=0.2) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  scale_color_manual(values = c("<2kbp" = "#F21A00", ">2kbp" = "#3B9AB2")) + 
  theme_classic() + facet_wrap(~cluster,scales="free",nrow = 1) + theme(axis.text = element_text(size = 20,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 15, angle = 0, hjust = 1, vjust = 0), axis.title.x = element_text(size = 20, angle = 0, hjust = .5, vjust = 0), axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),legend.text=element_text(size=20),legend.title=element_text(size=20))
dev.off()

#opn vs close
colnames(atac_broad)[169] <- "case.frq" 
colnames(atac_broad)[170] <- "control.frq" 
#atac_broad$cluster <- atac_broad$cluster_id
#atac_all <- rbind(atac_sub, atac_broad)
#atac_all$direction <- ifelse(atac_all$logFC_treatment > 0, 1, 0)
#ggvenn(list("Open"=atac_all$gene[atac_all$logFC_treatment > 0],"Closed"=atac_all$gene[atac_all$logFC_treatment < 0]),fill_color = c("#0073C2FF", "#F21A00"))

atac_broad$direction <- ifelse(atac_broad$logFC_treatment > 0, 1, 0)
atac_sub$direction <- ifelse(atac_sub$logFC_treatment > 0, 1, 0)
#ggvenn(list("Open"=atac_all$gene[atac_all$logFC_treatment > 0],"Closed"=atac_all$gene[atac_all$logFC_treatment < 0]),fill_color = c("#0073C2FF", "#F21A00"))
pie(table(atac_broad$direction),labels=c("Closed","Open"),col=c("darkslategrey","goldenrod2"),cex=2)
pdf("~/Desktop/atac_2023_paper/broad_openclose.pdf")
pie(table(atac_broad$direction),labels=c("Closed","Open"),col=c("darkslategrey","goldenrod2"),cex=2)
dev.off()

pie(table(atac_sub$direction),labels=c("Closed","Open"),col=c("darkslategrey","goldenrod2"),cex=2)
pdf("~/Desktop/atac_2023_paper/cluster_openclose.pdf")
pie(table(atac_sub$direction),labels=c("Closed","Open"),col=c("darkslategrey","goldenrod2"),cex=2)
dev.off()
##

#DAR volcano plots
#load peak-to-gene links
sub_p2g_full <- readRDS("~/Desktop/Subtype_p2g_redone.rds")
sub_p2g  <- sub_p2g_full %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45)
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE & atac$cluster_id=="ExN1",]
atac_sub$peakName <- atac_sub$gene
merge(atac_sub, sub_p2g, by="peakName") %>% data.frame () -> d
ggplot(d) +
  geom_point(aes(x = logFC_treatment, y = -log10(p_adj.loc_treatment), colour = sign(logFC_treatment))) + scale_colour_gradientn(colors = c("darkslategrey","goldenrod2")) + 
  geom_text_repel(aes(x = logFC_treatment, y = -log10(p_adj.loc_treatment), label = geneName),max.overlaps = 20) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") + theme_classic() + theme(legend.position = "none",plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = 20),axis.text=element_text(size=20))

atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE & atac$cluster_id=="Mic2",]
atac_sub$peakName <- atac_sub$gene
merge(atac_sub, sub_p2g, by="peakName") %>% data.frame () -> d
ggplot(d) +
  geom_point(aes(x = logFC_treatment, y = -log10(p_adj.loc_treatment), colour = sign(logFC_treatment))) + scale_colour_gradientn(colors = c("darkslategrey","goldenrod2")) + 
  geom_text_repel(aes(x = logFC_treatment, y = -log10(p_adj.loc_treatment), label = geneName),max.overlaps = 20) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") + theme_classic() + theme(legend.position = "none",plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = 20),axis.text=element_text(size=20))

#psygnet 
atac_subcluster <- atac[atac$logFC_treatment > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE,] 
sub_p2g$geneName[sub_p2g$peakName %in% atac_subcluster$gene] -> atac_genes
m1 <- psygenetGene(gene= unique(broad_genes),database = "ALL", verbose  = TRUE)
geneAttrPlot( m1, type = "evidence index" ) -> psy
saveRDS(psy,"~/Desktop/DownstreamResults/psygnet_plot.rds")
psy <- readRDS("~/Desktop/DownstreamResults/psygnet_plot.rds") #theme_ArchR(legendTextSize = 10) + 
ggplot(psy$data, aes(x = Category, y = value, fill = variable)) + geom_bar(stat = "identity", width = .6) + scale_fill_npg() + coord_flip() + ylab("Gene Disease Associations") + theme_classic() + theme(axis.text = element_text(size = 20))
#pdf("~/Desktop/atac_2023_paper/broad_psygnet_npg.pdf")
#ggplot(psy$data, aes(x = Category, y = value, fill = variable)) + geom_bar(stat = "identity", width = .6) + scale_fill_npg() + coord_flip() + ylab("Gene Disease Associations") + theme_classic() + theme(axis.text = element_text(size = 20))
dev.off()


