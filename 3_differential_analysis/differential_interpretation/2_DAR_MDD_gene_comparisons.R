
library(ggplot2)
library(psygenet2r)
library(tidyverse)
library(ggsci)

#load broad DARs
atac2 <- readRDS("~/res_freq_broad.rds")
tmp   <- atac2[, grep("Case.frq|Control.frq", colnames(atac2))]
tmp$Greater1pct <- ifelse(tmp$Case.frq  > 0.01 | tmp$Control.frq > 0.01, "TRUE", "FALSE")
atac2$Greater1pct <- tmp$Greater1pct
atac_broad <- atac2[abs(atac2$logFC_treatment) > log2(1.1) & atac2$p_adj.loc_treatment < 0.05 & atac2$Greater1pct==TRUE,] 
#load cluster DARs
atac <- readRDS("~/res_freq_subcluster.rds")
atac$cluster <- plyr::mapvalues(atac$cluster_id, from = unique(atac$cluster_id), to = c("Ast", "Ast", "Ast", "Ast", "Mic", "ExN", "ExN","ExN", "ExN","ExN", "ExN","ExN", "ExN", "ExN", "ExN", "InN", "InN", "InN", "InN", "Mic", "Mic", "Oli", "Oli", "Oli", "Oli", "Oli", "OPC", "OPC", "OPC"))
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater3pct==TRUE,]

#load linkages
sub_p2g_full <- readRDS("~/Subtype_p2g_redone.rds")
broad_p2g_full <- readRDS("~/broad_p2g.rds")
broad_p2g_full$peakName <- sub(':', '-', broad_p2g_full$peakName)
broad_p2g  <- broad_p2g_full %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45)
sub_p2g  <- sub_p2g_full %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45)

#DAR genes
sub_p2g$geneName[sub_p2g$peakName %in% atac_sub$gene] -> sub_genes
broad_p2g$geneName[broad_p2g$peakName %in% atac_broad$gene] -> broad_genes

####psygenet analysis####
m1 <- psygenetGene(gene= unique(c(broad_genes,sub_genes)),database = "ALL", verbose  = TRUE)
psy <- geneAttrPlot( m1, type = "evidence index" ) 
saveRDS(psy,"~/psygnet_plot.rds")
psy <- readRDS("~/psygnet_plot.rds")  
pdf("~/psygnet_plot_broad_sub.pdf",width=10)
ggplot(psy$data, aes(x = Category, y = value, fill = variable)) + geom_bar(stat = "identity", width = .6,position="dodge") + scale_fill_npg() + ylab("Gene Disease Associations") + theme_classic() + theme(axis.text = element_text(size = 20,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0), axis.title.x = element_text(size = 20, angle = 0, hjust = .5, vjust = 0), axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),legend.text=element_text(size=20),legend.title=element_blank(),legend.position="bottom") + ylim(0,180) 
dev.off()

##############################################################################################################
#GeneOverlap analysis

library(GenomicRanges)
library(plyranges)
library(GeneOverlap)

##compare DAR linked genes with snRNA
rna_subcluster <- read_csv("~/9_combination_pvals_all_unfilteredsubtype.csv")
#rna_subcluster$cluster <- plyr::mapvalues(rna_subcluster$cluster_id.male, from = unique(rna_subcluster$cluster_id.male), to = c("InN", "InN", "ExN_deep", "ExN_deep", "ExN_deep", "ExN_deep","ExN_uniform", "ExN_upper", "ExN_deep","ExN_uniform", "ExN_upper", "ExN_upper","ExN_upper", "ExN_deep",  "Ast", "Ast", "InN", "InN", "InN", "InN", "ExN_deep", "InN", "InN", "InN", "InN", "Oli", "Oli", "Oli", "OPC", "OPC", "ExN_uniform", "ExN_deep", "Mix", "End", "Mic"))
#rna_subcluster$cluster <- plyr::mapvalues(rna_subcluster$cluster_id.male, from = unique(rna_subcluster$cluster_id.male), to = c("InN", "InN", "ExN", "ExN", "ExN", "ExN","ExN", "ExN", "ExN","ExN", "ExN", "ExN","ExN", "ExN",  "Ast", "Ast", "InN", "InN", "InN", "InN", "ExN", "InN", "InN", "InN", "InN", "Oli", "Oli", "Oli", "OPC", "OPC", "ExN", "ExN", "Mix", "End", "Mic"))
rna_broad <- read_csv("~/9_combination_pvals_all_unfilteredbroad.csv")
##Filter DEGs at fdr 10% & logFC 
rna_subcluster$logFC<- (rna_subcluster$logFC.female + rna_subcluster$logFC.male) / 2
rna_broad$logFC<- (rna_broad$logFC.female + rna_broad$logFC.male) / 2
rna_subcluster_de <- rna_subcluster[rna_subcluster$signs > 0 & rna_subcluster$adjpval < 0.05 & abs(rna_subcluster$logFC) > log2(1.1) ,]
rna_broad_de <- rna_broad[rna_broad$signs > 0 & rna_broad$adjpval < 0.05 & abs(rna_broad$logFC) > log2(1.1) ,]
rna_deg <- c(rna_subcluster_de$gene, rna_broad_de$gene)

#MDD GWAS genes
#MDD Howard et al (2019)
gwas_genes <- read.csv("~/DownstreamResults/MDD_genes.csv")
gwas_howard <- c(gwas_genes$Howard2018, gwas_genes$Howard2019)
##gwas_howard <- gwas_howard[gwas_howard %ni% c("MDD_1","MDD_2","MDD_3", "")]

#MDD Als et al (2023)
gwas_als <- readxl::read_excel("~/Desktop/Als_magma.xlsx",col_names = "genes")

#psygenet
psy <- read.csv("~/Desktop/DownstreamResults/psygenet.xlsx.txt", sep="\t")
psy <- psy[psy$PsychiatricDisorder %in% c("Depressive disorders"),]

##run GeneOverlap analysis
Genes <- list(GWAS.Howard=gwas_howard, GWAS.Als=gwas_als$genes, PsyGeNet=psy$Gene_Symbol,snRNA_DEGs=rna_deg) 
ATAC <- list(broad_DAR=broad_genes, subcluster_DAR=sub_genes) 
go.obj <- newGOM(rna_list, ATAC,genome.size=length(unique(gc$nearestGene)))

##Plot
pdf("~/mdd_atac_genes_comparison.pdf")
drawHeatmap(go.obj, adj.p = T,log.scale = T,ncolused = 9)
dev.off()


