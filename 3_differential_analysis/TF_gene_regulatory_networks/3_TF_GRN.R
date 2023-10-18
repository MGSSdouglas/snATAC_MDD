##draw cell-type TF-gene network for DARs##
##@inspired by morabito et. al (2021) TF-gene network##
##https://github.com/swaruplab/Single-nuclei-epigenomic-and-transcriptomic-landscape-in-Alzheimer-disease/blob/master/snATAC_TF-networks.Rmd

library(Matrix)
library(GenomicRanges)
library(tidyverse)
library(ggplot2)
library(igraph)
library(qgraph)
library(ArchR)
library(GenomicRanges)
library(GeneOverlap)
library(clusterProfiler)
library(org.Hs.eg.db)

##load all subcluster peaks
peaks <- readRDS("~/Desktop/IterativelyRemovedPeaks_sub.rds")
peaks$peaks <- paste(seqnames(peaks),ranges(peaks),sep="-")

#snRNA-seq combined model approach
#rna_subcluster <- read.csv("~/Desktop/Combined_model_unfilteredCluster.csv")
#rna_broad <- read.csv("~/Desktop/Combined_model_unfilteredBroad.csv")
#rna_subcluster$id <- plyr::mapvalues(rna_subcluster$cluster_id, from = unique(rna_subcluster$cluster_id), to = c("InN", "InN", "ExN_deep", "ExN_deep", "ExN_deep", "ExN_deep","ExN_uniform", "ExN_upper", "ExN_deep","ExN_uniform", "ExN_upper", "ExN_upper","ExN_upper", "ExN_deep",  "Ast", "Ast", "InN", "InN", "InN", "InN", "ExN_deep", "InN", "InN", "InN", "InN", "Oli", "Oli", "Oli", "OPC", "OPC", "ExN_uniform", "ExN_deep", "Mix", "End", "Mic"))
#rna_subcluster_de <- rna_subcluster[rna_subcluster$p_adj.loc < 0.1 & abs(rna_subcluster$logFC) > log2(1.1) & rna_subcluster$Greater3==TRUE,]
#rna_broad_de <- rna_broad[rna_broad$p_adj.loc < 0.1 & abs(rna_broad$logFC) > log2(1.1) & rna_broad$Greater3==TRUE,]

##snRNA-seq meta-analysis with p_val combination approach
rna_broad <- read_csv("~/Desktop/9_combination_pvals_all_unfilteredbroad.csv")
rna_subcluster <- read_csv("~/Desktop/9_combination_pvals_all_unfilteredsubtype.csv")
rna_subcluster$id <- plyr::mapvalues(rna_subcluster$cluster_id.male, from = unique(rna_subcluster$cluster_id.male), to = c("InN", "InN", "ExN_deep", "ExN_deep", "ExN_deep", "ExN_deep","ExN_uniform", "ExN_upper", "ExN_deep","ExN_uniform", "ExN_upper", "ExN_upper","ExN_upper", "ExN_deep",  "Ast", "Ast", "InN", "InN", "InN", "InN", "ExN_deep", "InN", "InN", "InN", "InN", "Oli", "Oli", "Oli", "OPC", "OPC", "ExN_uniform", "ExN_deep", "Mix", "End", "Mic"))
rna_subcluster$cluster <- plyr::mapvalues(rna_subcluster$cluster_id.male, from = unique(rna_subcluster$cluster_id.male), to = c("InN", "InN", "ExN", "ExN", "ExN", "ExN","ExN", "ExN", "ExN","ExN", "ExN", "ExN","ExN", "ExN",  "Ast", "Ast", "InN", "InN", "InN", "InN", "ExN", "InN", "InN", "InN", "InN", "Oli", "Oli", "Oli", "OPC", "OPC", "ExN", "ExN", "Mix", "End", "Mic"))

##Filter DEGs at fdr 10% & logFC & signs==1 (same direction effects in males and females)
rna_subcluster$logFC<- (rna_subcluster$logFC.female + rna_subcluster$logFC.male) / 2
rna_broad$logFC<- (rna_broad$logFC.female + rna_broad$logFC.male) / 2
rna_subcluster_de <- rna_subcluster[rna_subcluster$signs > 0 & rna_subcluster$adjpval < 0.1 & abs(rna_subcluster$logFC) > log2(1.1) ,]
rna_broad_de <- rna_broad[rna_broad$signs > 0 & rna_broad$adjpval < 0.1 & abs(rna_broad$logFC) > log2(1.1) ,]
rna_de_up <- c(rna_subcluster_de$gene[rna_subcluster_de$logFC > 0], rna_broad_de$gene[rna_broad_de$logFC > 0])
rna_de_down <- c(rna_subcluster_de$gene[rna_subcluster_de$logFC < 0], rna_broad_de$gene[rna_broad_de$logFC < 0])

#Gandal et. al #bulk MDD DEGs
gandal_genes <- read.csv("~/Desktop/Gandal_MDD_genes.csv",header=T)
gandal_genes <- unique(gandal_genes$external_gene_id[gandal_genes$MDD.FDR < 0.1])

##list of GWAS, HiC GWAS genes (https://github.com/thewonlab/H-MAGMA) and psygenet genes 
gwas_genes <- read.csv("~/Desktop/DownstreamResults/MDD_genes.csv")
gwas_mdd_wray_howard <- c(gwas_genes$Howard2018, gwas_genes$Howard2019, gwas_genes$Wray2018)
gwas_mdd_wray_howard <- gwas_mdd_genes[gwas_mdd_genes %ni% c("MDD_1","MDD_2","MDD_3", "")]
gwas_mdd_silveria <- readxl::read_excel("~/Desktop/Silveria_genes.xlsx",col_names = "genes")
gwas_mdd_als <- readxl::read_excel("~/Desktop/Als_magma.xlsx",col_names = "genes")
gwas_mdd_levey <- c("ADARB2","CCDC71","DDX27","GPR27","HISTIH2BN","PARK2","SELENOH","SLC25A17","SORBS3","ST13","ZKSCAN4","ZNF502","ACVRL1","CKB","CTTNBP2","L3MBTL2","MSH5","PCDHA8","RCN1","RP11-1I2.1","RP11-24H2.3","RP11-3B7.1","RP1-86C11.7","SF3B1","SMIM4","WDR55","CELF4","CHADL","FCF1","LINC01415","LINC02060","PAX6","RAB27B","RASGRP1","ZNF184","BTN2A1","CTA-14H9.5","CTC-467M3.3","EBLN3P","FADS1","HCG11","HSPA1A","PCDHA1","PCDHA2","PCDHA3","PGBD1","POU5F2","RBX1","TAL1","UBA7","ZNF852","AREL1","ASCC3","BSN","DENND1B","FADS2","KLC1","KLF7","OTOA","PRSS16","STK24","TMEM42","ZKSCAN7","ZNF501","ZSCAN9","BAG5","CRB1","CTD-2298J14.2","HSPA1L","LIN28B","LRFN5","PROX2","RP1-265C24.8","STAU1","TCAIM","TMEM161B","YLPM1","APEH","CENPL","CLOCK","FAM120A","PCDHA5","PPP3CC","SPPL3","TRMT61A","ZDHHC5","B4GALNT4","BSN-AS2","LINC00461","NT5C2","RHOA","RP1-153G14.4","SLC12A5","TIAF1","ZSCAN12","CSE1L","EP300","FURIN","ITPR3","MLF1","NEGR1","NICN1","NLGN1","RP11-220I1.5","TRAF3","ZKSCAN8","BTN3A3","CCDC36","CD40","DRD2","SORCS3","TMEM106B","USP15","USP4","ZC3H7B","ZNF197","ZNF445","AMT","FAM120AOS","GPX1","LIN28B-AS1","TRIM8","HIST1H1B","HIST1H2BO","KIAA1143","KLHDC8B","LINC00997","PCDHA7","RP11-571M6.17","TMEM258","BTN3A1","C15orf59-AS1","CTC-498M16.4","DAG1","EIF4E3","RP11-318C24.2","ZNFX1","ZSCAN16","ANKK1","BTN2A2","BTN3A2","CITF22-92A6.1","CLP1","DCP2","KCNG2","LAMB2","PTCH1","RERE","SCRN3","TCTA","TP53I11","ZNF165","ZNF660","SORCS3","RP11-467L20.10","DAGLA","FADS2","TMEM258","MYRF","FADS3","RAB3IL1","FEN1","FADS1","MIR611","MIR1908","NR4A1","ACVR1B","GRASP","ANKRD33","RNU6-574P","ACVRL1","RP11-1100L3.4","ACVR1B","GRASP","ANKRD33","RNU6-574P","ACVRL1","RP11-1100L3.4","HNF1A","C12orf43","SPPL3","ARF1P2","CLIC1P1","HNF1A-AS1","RPL12P33","AC079602.1","RP11-216P16.2","AL450423.1","PCDH8P1","RN7SL618P","ESR2","SYNE2","YLPM1","AC007956.1","FCF1","AREL1","SNORA7","RP11-173A8.2","AMN","TRAF3","CDC42BPB","RNU6-1316P","AL117209.1","snoU13","RP11-661D19.3","RP11-365N19.2","PHF12","SEZ6","PIPOX","MYO18A","TIAF1","AC024619.2","RP11-321A17.3","RP11-321A17.5","RP11-321A17.4","PHF12","SEZ6","PIPOX","MYO18A","TIAF1","AC024619.2","RP11-321A17.3","RP11-321A17.5","RP11-321A17.4","CELF4","MIR4318","AC011225.1","RP11-188I24.1","RP11-142I20.1","VRK2","CTD-2026C7.1","CCDC36","LAMB2","USP19","QARS","CCDC71","KLHDC8B","C3orf62","QRICH1","C3orf84","RP11-3B7.1","RN7SL182P","RP11-694I15.7","RP11-3B7.7","ARAP2","RP11-6N13.1","RP11-284A20.1","RP11-284A20.2","RP11-284A20.3")
mdd_hmagma <- readxl::read_excel("~/Desktop/hmagma_adult_mdd.xlsx")
mdd_hmagma$fdr <- p.adjust(mdd_hmagma$P, method = "BH")
library('biomaRt') #convert genes from ensembl id to symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "hgnc_symbol"),values=mdd_hmagma$GENE,mart= mart)
mdd_hmagma_genes <- mdd_hmagma$GENE[mdd_hmagma$fdr < 0.1]
mdd_hmagma_genes <- G_list$hgnc_symbol[G_list$ensembl_gene_id %in% mdd_hmagma_genes]

#all gwas mdd genes
gwas_mdd_genes <- c(gwas_mdd_wray_howard, levey, als$genes, silveria$genes, mdd_hmagma_genes)

psy <- read.csv("~/Desktop/DownstreamResults/psygenet.xlsx.txt", sep="\t")
psy <- psy[psy$PsychiatricDisorder %in% c("Depressive disorders"),]

##significant p2g linkages
sub_p2g_full <- readRDS("~/Desktop/Subtype_p2g_redone.rds")
sub_p2g  <- sub_p2g_full %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45)


cluster=c("Mic") #
cluster=c("ExN_deep")
cluster="ExN"
cluster="Ast" #, "ExN_uniform", "ExN_upper")

##Atac DARs filter
atac <- readRDS("~/Desktop/res_freq_subcluster.rds")
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
#tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.05 | tmp$control.frq > 0.05, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct
#use 0.01 fdr for mic network as too large
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$cluster_id %in% "ExN1" & atac$Greater3pct==TRUE,] #| abs(atac$logFC_male.trt) > log2(1.1) & atac$p_adj.loc_male.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_female.trt) > log2(1.1) & atac$p_adj.loc_female.trt < 0.2 & atac$Greater3pct==TRUE | abs(atac$logFC_interaction) > log2(1.1) & atac$p_adj.loc_interaction < 0.2 & atac$Greater3pct==TRUE,]
atac_top <- atac_sub[order(atac_sub$p_adj.loc_treatment),]
atac_top$peakName <- atac_top$gene
join <- left_join(atac_top, sub_p2g, by="peakName")
top_atac_genes <- join %>% group_by(peakName) %>% filter(Correlation==max(Correlation)) #plot some top most correlated genes
top_atac_genes <- c(top_atac_genes$geneName[1:20])

#top_atac_genes <- top_atac_genes[top_atac_genes %ni% c("LRRC8D","SLC35D2","CPED1")] #for mic; remove some as were overlapping
#use 0.01 fdr for mic network as too large


##load motif data for homer
#motif_names <- readRDS("~/Desktop/DownstreamResults/JASPAR_motif_names.rds")
motif_matrix <- readRDS("~/Desktop/motif_matrix_homer.rds")
motif_positions <- readRDS("~/Desktop/motif_positions_homer.rds")


motifs <- c("Bach2.bZIP_17", "CTCF.Zf_43","Ascl1.bHLH_9","X.box.HTH_304") #for ExN1
motifs <- c("PU.1.ETS_235","IRF1.IRF_136","Atf7.bZIP_14","Sp1.Zf_266")   #for Mic2
#motifs <- c("GRE.NR..IR3_117","GRE.NR..IR3_118", "AR.halfsite.NR_7","Bach2.bZIP_17","CTCF.Zf_43")
type = c("both") 
peaks_vector <- vector()
edge_df <- data.frame()
vertex_df <- data.frame()
cur_motif <- vector()
#shared <- peaks_vector[duplicated(peaks_vector)]  
for(motif in motifs){
  print(motif)
  cur <- do.call(rbind,str_split(motif, "_"))[1]
  motifs_tmp <- motif_matrix[,colnames(motif_matrix) %in% motif]
  motifs_tmp <- motifs_tmp[motifs_tmp >= 1]
  tmp <- peaks[peaks$peaks %in% names(motifs_tmp),]
  
  #subset to DAR peaks 
  tmp2 <- tmp[tmp$peaks %in% c(atac_sub$gene),] 
  
  #proximal genes
  nearest_genes <- c(sub_p2g$geneName[sub_p2g$peakName %in% tmp2$peaks & sub_p2g$peakType %in% c("Promoter")],tmp2$nearestGene[tmp2$peakType %in% c("Promoter")])
  #print(length(nearest_genes))
  #distal genes
  enhancer_genes <- unique(c(sub_p2g$geneName[sub_p2g$peakName %in% tmp2$peaks & sub_p2g$peakType %in% c("Distal","Intronic")])) #,tmp2$nearestGene[tmp2$peakType %in% c("Distal","Intronic")])) #tmp2$peaks & sub_p2g$peakType == "Distal"])
  
  #print TF DAR peaks and genes 
  print(length(tmp2))
  print(length(unique((c(nearest_genes,enhancer_genes)))))
  
  if(length(nearest_genes) > 0){
    cur_promoter_edge_df <- data.frame(
      from=cur,
      to=as.character(unlist(nearest_genes)),
      type='Promoter'
    )
  } else{cur_promoter_edge_df <- data.frame()}
  
  if(length(enhancer_genes) > 0){
    cur_enhancer_edge_df <- data.frame(
      from=cur,
      to=as.character(unlist(enhancer_genes)),
      type='Distal'
    )
  } else{cur_enhancer_edge_df <- data.frame()}
  #edge_df <- rbind(cur_promoter_edge_df[1:10], cur_distal_edge_df[1:10])
  if(type=="Promoter"){
    cur_vertex_df <- data.frame(name = c(cur, as.character(unique(c(unlist(nearest_genes))))))
    cur_edge_df <- rbind(cur_promoter_edge_df)
    #cur_edge_df$from <- c(cur,length(unique(c(nearest_genes,enhancer_genes))))
    
  } else if(type=="Distal"){
    cur_vertex_df <- data.frame(name = c(cur,as.character(unique(c(unlist(enhancer_genes))))))
    cur_edge_df <- rbind(cur_enhancer_edge_df)
    #cur_edge_df$from <- c(paste0(cur, " (",length(unique(tmp2$peaks)),", ",length(unique(c(nearest_genes,enhancer_genes))),")"))
    
  }else{
    cur_vertex_df <- data.frame(name = c(cur,as.character(unique(c(unlist(nearest_genes), unlist(enhancer_genes))))))
    cur_edge_df <- rbind(cur_promoter_edge_df, cur_enhancer_edge_df)
    cur_edge_df$from <- c(paste0(cur))
    
  }
  
  edge_df <- rbind(edge_df,cur_edge_df)
  vertex_df <- rbind(vertex_df, cur_vertex_df)
  #cur <- c(paste0(cur)) #" (",length(unique(tmp2$peaks)),", ",length(unique(c(nearest_genes,enhancer_genes))),")"))
  cur_motif <- c(cur_motif, cur)
  #overlapping peaks
  #peaks_vector <- c(peaks_vector, tmp2$peaks) 
}


vertex_df <- data.frame(name=na.omit(as.character(unique(vertex_df$name))))
vertex_df$name <- as.character(vertex_df$name)
edge_df <- na.omit(edge_df)

#shared_genes <- c(shared_nearest_genes,shared_enhancer_genes)
#edge_df$type <- ifelse(edge_df$to %in% shared_genes, "shared", edge_df$type)

g <- igraph::graph_from_data_frame(edge_df, directed=TRUE, vertices=vertex_df)


################################################################################
# visual settings for network
################################################################################


# remove labels if gene is not DE, or not a TF:
up <-  unique(rna_de_up)  
down <- unique(rna_de_down) 
de_targets <- as.character(vertex_df$name[vertex_df$name %in% unique(c(up,down))])
all_genes <- c(gwas_mdd_genes,psy$Gene_Symbol,gandal_genes,de_targets)
mdd_genes <- c(de_targets,gandal_genes)
gwas_genes <- c(gwas_mdd_genes,psy$Gene_Symbol)

vertex_df$label <- ifelse(vertex_df$name %in%  top_atac_genes, vertex_df$name,"  ") #vertex_df$label for ast OR ""
vertex_df$label <- ifelse(vertex_df$name %in% c(as.character(unique(all_genes))), vertex_df$name, vertex_df$label)
vertex_df$label <- ifelse(vertex_df$name %in% c(as.character(motif_names),cur_motif), vertex_df$name,  vertex_df$label)

# set node color based on DEGs:
vertex_df$color <- ifelse(vertex_df$name %in% c(as.character(motif_names),cur_motif), 'darkseagreen', 'ghostwhite')
#vertex_df$color <- ifelse(vertex_df$name %in% up, "gold2", vertex_df$color)
vertex_df$color <- ifelse(vertex_df$name %in% mdd_genes, 'darkorange', vertex_df$color) #cadetblue1
vertex_df$color <- ifelse(vertex_df$name %in% gwas_genes, 'hotpink', vertex_df$color)
vertex_df$color <- ifelse(vertex_df$name %in% intersect(gwas_genes,mdd_genes), 'cornflowerblue', vertex_df$color)
# italics font for genes:
vertex_df$font <- ifelse(vertex_df$name %in% as.character(motif_names), 2, 4)

# set size to larger if the gene is a TF:
vertex_df$size <- ifelse(vertex_df$name %in% c(as.character(motif_names),cur_motif), 8, 0.2)
vertex_df$size <- ifelse(vertex_df$name %in% all_genes, 5, vertex_df$size)


promoter_color <- 'goldenrod2'
enhancer_color <-  'darkcyan'         #'lavenderblush2' #'mistyrose2' #
shared_color <- 'powderblue'


edge_colors <- ifelse(E(g)$type == 'Promoter', promoter_color,ifelse(E(g)$type == 'Distal', enhancer_color, shared_color))
V(g)$label.cex <- 0.6

set.seed(4) #38 for ExN #3 for mic  #27, 004 for ast 
l <- layout_with_fr(g)
#svg("~/Desktop/atac_2023_paper/TF_network_selection/ExN1_network.svg", width=10, height=8)
pdf("~/Desktop/atac_2023_paper/TF_final_selection/ExN1_network.pdf")
plot(
  g, layout=l,
  vertex.size=vertex_df$size,
  edge.color=edge_colors,
  edge.alpha=0.5,vertex.shape=vertex_df$shape,
  vertex.color=vertex_df$color,vertex.frame.color='white',vertex.label.font=.000001,vertex.frame.width=0.1,
  vertex.label=vertex_df$label, vertex.label.family="Times", vertex.label.font=vertex_df$font,
  vertex.label.color = 'black',vertex.label.degree = pi/4,vertex.label.dist=0.1,
  edge.arrow.size=0.3
)
dev.off()


