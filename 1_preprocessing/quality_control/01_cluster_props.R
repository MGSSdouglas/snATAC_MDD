library(ArchR)
library(ggplot2)
library(tidyverse)

#Load ArchR project
proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered")

#calculate cluster props
Subject_props <- table(proj$SubClusters, proj$Subject) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Subject_props", "Freq"))
Batch_props <- table(proj$SubClusters,  proj$batch) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Batch_props", "Freq"))
Sex_props <- table(proj$SubClusters, proj$sex) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Sex_props", "Freq"))
Condition_props <- table(proj$SubClusters, proj$condition) %>% as.data.frame() %>% setNames(nm = c("Cluster", "Condition_props", "Freq"))

#plot
pdf("~/scratch/snATAC_plots_to_submit/subject_props.pdf",height=10, width=20)
ggplot(Subject_props, aes(fill=Subject_props, y=Freq, x=Cluster)) + geom_bar(position="fill", stat="identity")+ 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20),axis.text.y = element_text(size = 20))
dev.off()

pdf("~/scratch/snATAC_plots_to_submit/sex_props_corrected.pdf",height=3, width=10) 
ggplot(Sex_props, aes(fill=Sex_props, y=Freq, x=Cluster)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=c("brown","cornflowerblue")) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),axis.text.y = element_text(size = 10))
dev.off()

proj$Age <- cut(proj$Age,breaks=c(0, 20, 40, 60, 80, 100),labels=c('0-20', '20-40', '40-60', '60-80', '80-100'))
age_props <- table(proj$SubClusters,  proj$Age) %>% as.data.frame() %>% setNames(nm = c("Cluster", "age_props", "Freq"))
ggplot(age_props, aes(fill=age_props, y=Freq, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20),axis.text.y = element_text(size = 20))

proj$PMI <- cut(proj$PMI,breaks=c(0, 20, 40, 60, 80, 100),labels=c('0-20', '20-40', '40-60', '60-80', '80-100'))
PMI_props <- table(proj$SubClusters,  proj$PMI) %>% as.data.frame() %>% setNames(nm = c("Cluster", "PMI_props", "Freq"))
pdf("~/scratch/snATAC_plots_to_submit/PMI_props_corrected.pdf",width=12)
ggplot(PMI_props, aes(fill=PMI_props, y=Freq, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20),axis.text.y = element_text(size = 20))
dev.off()

my_theme <- theme(	
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.position="bottom",
  plot.title=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)
pdf("~/scratch/snATAC_plots_to_submit/plot_QC/umap_age.pdf")
p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Age", embedding = "UMAPHarmony",baseSize=20,labelMeans=FALSE, legendSize =7, labelSize = 10,labelAsFactors = FALSE) + my_theme 
dev.off()
pdf("~/scratch/snATAC_plots_to_submit/plot_QC/umap_sex.pdf")
p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sex", embedding = "UMAPHarmony",baseSize=20,labelMeans=FALSE, legendSize =7, labelSize = 10,labelAsFactors = FALSE) + my_theme 
dev.off()

