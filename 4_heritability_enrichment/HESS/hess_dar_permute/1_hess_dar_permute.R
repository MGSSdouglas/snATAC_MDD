library(GenomicRanges)
library(plyranges)
library(Signac)
library(tidyverse)
library(biomaRt)
library(ArchR)
library(reshape2)

#load hess LD blocks
hess <- read.csv("~/Desktop/step2.txt",se="\t")
hess$peaks <- paste(hess$chr, hess$start, hess$end, sep="-")
hess$peaks <- paste0("chr",hess$peaks)
hess <- hess[order(hess$local_h2g, decreasing = TRUE),]
hess_gr <- Signac::StringToGRanges(hess$peaks)
#hess2 <- hess[hess$local_h2g > 0,]
hess_sig <- hess[hess$p < 0.05,] #90 LD blocks with significant MDD heritability

#differential candidate peaks in MDD
atac <- readRDS("~/Desktop/res_freq_subcluster.rds")
atac$cluster <- plyr::mapvalues(atac$cluster_id, from = unique(atac$cluster_id), to = c("Ast", "Ast", "Ast", "Ast", "Mic", "ExN", "ExN","ExN", "ExN","ExN", "ExN","ExN", "ExN", "ExN", "ExN", "InN", "InN", "InN", "InN", "Mic", "Mic", "Oli", "Oli", "Oli", "Oli", "Oli", "OPC", "OPC", "OPC"))
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater3pct <- ifelse(tmp$case.frq  > 0.03 | tmp$control.frq > 0.03, "TRUE", "FALSE")
atac$Greater3pct <- tmp$Greater3pct
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.2 & atac$Greater3pct==TRUE,] 

#all open chromatin peaks
allPeaks <- readRDS("~/Desktop/IterativelyRemovedPeaks_sub.rds")

#differential peaks enriched in LD blocks of high heritability than random peaks in any LD block 
#loop over every cluster
clusters <- unique(atac_sub$cluster_id)
test.id <- rep(0, length(clusters))
for(cluster in clusters){
  print(cluster)
  tmp_cluster <- atac_sub[atac_sub$cluster_id==cluster,]
  #no. of DAR peaks within significant LD blocks
  test.id <- length(subsetByOverlaps(Signac::StringToGRanges(tmp_cluster$gene), Signac::StringToGRanges(hess_sig$peaks), type="within"))/length(tmp_cluster$gene)
  assign(paste0("test",cluster),test.id)
  
  #set-up permutations for peaks 
  perm = 10000 #no. of permutations
  P <- perm
  n_peaks <- length(tmp_cluster$gene) #no. of peaks in that cluster
  
  #retain peaks in same chromosomes as test peaks.
  variable_gr <- allPeaks %>% filter(seqnames %in% levels(seqnames(Signac::StringToGRanges(tmp_cluster$gene))))
  variable_gr$peaks <- paste(seqnames(variable_gr),ranges(variable_gr),sep="-")
  PermSamples <- matrix(0, nrow = n_peaks, ncol = P)
  for(i in 1:P)
  {
    PermSamples[, i] <- sample(variable_gr$peaks, 
                               size = n_peaks, 
                               replace = FALSE)
  }
  #any LD block
  n_ld <- length(hess_sig$peaks)
  variable_ld <- hess$peaks
  PermSamples2 <- matrix(0, nrow = n_ld, ncol = P)
  for(i in 1:P)
  {
    PermSamples2[, i] <- sample(variable_ld, 
                                size = n_ld, 
                                replace = FALSE)
  }
  Perm.peaks <-  rep(0, P)
  for (i in 1:P)
  {
    # no. of sampled peaks within any LD block
    print(i)
    Perm.peaks[i] <- length(subsetByOverlaps(Signac::StringToGRanges(PermSamples[,i]),Signac::StringToGRanges(PermSamples2[, i]),type="within"))/length(tmp_cluster$gene)
    assign(paste0("perm",cluster),Perm.peaks)
  }
}
#print emperical significance
for(cluster in unique(atac_sub$cluster_id)){
  print(cluster)
  print(mean(get(paste0("perm",cluster)) >= get(paste0("test",cluster))))
}

#plot
#pie(c(testAst1,testAst2,testAst3,testAst4,testEnd1,testExN1,testExN1_L23,testExN1_L24,testExN1_L56,testExN2_L23,testExN3_L56,testMic1,testMic2,testOli2,testOli4,testOli5,testOli6,testOli7,testOPC1), labels=c("Ast1","Ast2","Ast3","Ast4","Mic3","ExN1_L36","ExN1_L23","ExN1_L24","ExN1_L56","ExN2_L23","ExN3_L56","Mic1","Mic2","Oli2","Oli4","Oli5","Oli6","Oli7","OPC1"))
test <- c(testAst1,testAst2,testAst3,testAst4,testExN1, testExN1_L23, testExN1_L24, testExN1_L56, testExN2, testExN2_L23, testExN2_L46, testExN2_L56, testExN3_L46, testExN3_L56, testIn_LAMP5,testIn_PV, testIn_VIP, testMic1, testMic2,testEnd1, testOli2, testOli4, testOli5, testOli6, testOli7, testOPC1, testOPC3)
df <- data.frame(permAst1,permAst2,permAst3,permAst4,permExN1, permExN1_L23, permExN1_L24, permExN1_L56, permExN2, permExN2_L23, permExN2_L46, permExN2_L56, permExN3_L46, permExN3_L56, permIn_LAMP5,permIn_PV, permIn_VIP, permMic1, permMic2,permEnd1, permOli2, permOli4, permOli5, permOli6, permOli7, permOPC1, permOPC3)
colnames(df) <- c("Ast1", "Ast2", "Ast3", "Ast4","ExN1", "ExN1_L23", "ExN1_L24", "ExN1_L56", "ExN2", "ExN2_L23", "ExN2_L46", "ExN2_L56", "ExN3_L46","ExN3_L56", "In_LAMP5", "In_PV", "In_VIP", "Mic1", "Mic2",  "Mic3","Oli2", "Oli4", "Oli5", "Oli6", "Oli7", "OPC1", "OPC3")
colnames(df) <- plyr::mapvalues(colnames(df), from = colnames(df), c("Ast1", "Ast2", "Ast3" , "Ast4" , "ExN1", "ExN11", "ExN10", "ExN3", "ExN2", "ExN12","ExN8", "ExN4", "ExN9",  "ExN5", "InLAMP5", "InPV", "InVIP", "Mic1", "Mic2", "End1", "Oli2", "Oli4", "Oli5", "Oli6", "Oli7", "OPC1", "OPC3"))

df2 <- melt(df)
df2$value <- df2$value *100
colnames(df2) <- c("Clusters", "Percentage")
test<-rep(test,each=10000)
test <- test * 100
pdf("~/Updated_boxplot_Hess_1k_permuted_DARs.pdf", width=10, height=6)
ggplot(df2, aes(x=Clusters,y=Percentage, fill=Clusters)) + geom_boxplot() + geom_point(data = df2,aes(x=Clusters, y=test),color = 'brown',position=position_dodge(0.8),shape=10, size=5) +  theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text = element_text(size=20),axis.title = element_text(size=20)) + ylim(0,30) + ylab("% DARs in MDD heritable LD blocks") + scale_fill_manual(values = c("#23B205","#74c69d","#2DF304","#145F04","#E7B321","#a6a3a1","#ccc7c4","#fab866","#FAE605","#736f6c", "#F7F7B5","#f29f3a","#E7D973","#d18221","#53B1F3","#4853A2","#6E7EBD","#81B2F1","#1792AA","#83D6E6","#E3D2A7","#B199C7","#7F56A4","#E04EC0","#E677CD","#F4BBE7", "#f4978e","#F8651B", "#f08080"))
dev.off()

