library(Seurat)
library(ggpubr)
library(tidyverse)
library(GeneOverlap)


#GM vs WM microglia: https://www.nature.com/articles/s41467-019-08976-7
#depDAM microglia MDD: https://www.biologicalpsychiatryjournal.com/article/S0006-3223(23)01239-8/fulltext

mic <- readRDS("~/scratch/cluster_profile/Mic_snRNA_GeneIntegrationMatrix_ArchR.rds")

Idents(mic) <- "Clusters"
mic <- subset(mic, idents=c("Mic1","Mic2"))

GM <- read.csv("~/scratch/cluster_profile/GM_WM_Microglia.csv")
WM <- read.csv("~/scratch/cluster_profile/WM_GM_Microglia.csv")
GM <- GM[order(GM$FoldChange,decreasing=T),]
WM <- WM[order(WM$FoldChange,decreasing=T),]

mic <- AddModuleScore(mic, features=list(GM=GM$X[1:25],WM=WM$X[1:25]))

#top 25 sorted by logFC 
pdf("~/scratch/cluster_profile/Mic_GM.pdf",height=5,width=5)
mic@meta.data %>%  ggplot() + aes(y = Cluster1, x = Clusters, fill= Clusters) + geom_boxplot() + stat_compare_means(method = "wilcox.test") + theme_classic() + theme(axis.text = element_text(size = 20),axis.title =element_text(size = 20)) + ylab("GM Module Score")  + scale_fill_manual(values=c("#E3D2A7","#1792AA","#83D6E6"))
dev.off()
pdf("~/scratch/cluster_profile/Mic_WM.pdf",height=5,width=5)
mic@meta.data %>%  ggplot() + aes(y = Cluster2 , x = Clusters, fill= Clusters) + geom_boxplot()  + stat_compare_means(method = "wilcox.test") + theme_classic() + theme(axis.text = element_text(size = 20),axis.title =element_text(size = 20)) + ylab("WM Module Score") + scale_fill_manual(values=c("#1792AA","#83D6E6")) #"#E3D2A7",
dev.off()

#MDD DAM #https://www.medrxiv.org/content/10.1101/2023.01.11.23284393v1.full#F7
dam <- read.csv("/home/anjali5/scratch/cluster_profile/MDD_DAM.csv")
dim(x)
#[1] 81  1
go.obj <- newGeneOverlap(dam$GENE,mic_genes,genome.size=length(unique(gc$nearestGene))) #gc has genes detected in atac ArchR
go.obj <- testGeneOverlap(go.obj)
print(go.obj)


mdd_down <- read.csv("/home/anjali5/scratch/cluster_profile/MDD_DAM.csv")
mdd_up <- c("LTBP3", "BCL9L", "C7orf26", "TFCP2L1", "SHROOM1", "LINCO1091", "CNTNAP2", "LIMCH1", "SCHIP1", "RP11-580I16.2", "PTPRB")
mic <- AddModuleScore(mic, features=list("MDD1"= mdd$GENE, "MDD2"= mdd_up), name="MDD",search=T)
mic@meta.data %>%  ggplot() + aes(y = MDD1 , x = Clusters, fill= condition) + geom_violin() + stat_compare_means() + theme_classic() +  theme(axis.text = element_text(size = 20),axis.title =element_text(size = 20))

pdf("~/scratch/cluster_profile/Down_DAM_Mic.pdf",height=4,width=4)
mic@meta.data %>%  ggplot() + aes(y = MDD1, x = Clusters, fill=SubClusters) + geom_boxplot() + stat_compare_means(method="wilcox.test") +  ylab("MDD DAM Downregulated") + theme_classic() +  theme(axis.text = element_text(size = 20),axis.title =element_text(size = 20)) + scale_fill_manual(values=c("#1792AA","#83D6E6"))
dev.off()
pdf("~/scratch/cluster_profile/Up_DAM_Mic.pdf",height=4,width=4)
mic@meta.data %>%  ggplot() + aes(y = MDD2, x = SubClusters, fill=SubClusters) + geom_boxplot() + stat_compare_means(method="wilcox.test") + ylab("MDD DAM Upregulated") +  theme_classic() +  theme(axis.text = element_text(size = 20),axis.title =element_text(size = 20)) + scale_fill_manual(values=c("#1792AA","#83D6E6"))
dev.off()


