library(GenomicRanges)
library(plyranges)
library(GeneOverlap)
library(plotrix)
library(ggplot2)
#####################################################################################
#load psychencode NeuN+/- ACC & PFC consolidated H3K27ac and H3k4me3 histone peaks 
#https://pubmed.ncbi.nlm.nih.gov/30038276/
#converted to hg38
promoter <- read_bed("/Users/anjalichawla/Desktop/Diff_histone_peaks/hg38_ncs/H3K4me3.bed")
enhancer <- read_bed("/Users/anjalichawla/Desktop/Diff_histone_peaks/hg38_ncs/H3K27ac.bed")
enhancer$peaks <- paste(seqnames(enhancer), ranges(enhancer), sep="-")
promoter$peaks <- paste(seqnames(promoter), ranges(promoter), sep="-")

#overlap with DARs open 
atac_up <- atac_all[atac_all$logFC_treatment > 0,]
subsetByOverlaps(Signac::StringToGRanges(c(atac_up$gene)), unique(c(enhancer))) -> x
x$peaks <- paste(seqnames(x),ranges(x),sep="-")
subsetByOverlaps(Signac::StringToGRanges(c(atac_up$gene)), unique(c(promoter))) -> y
y$peaks <- paste(seqnames(y),ranges(y),sep="-")
atac_up$histone <- ifelse(atac_up$gene %in% x$peaks, "H3K27ac", ifelse(atac_up$gene %in% y$peaks, "H3K4me3", ifelse(atac_up$gene %in% intersect(x$peaks, y$peaks), "H3K4me3+H3K27ac", "other")))
atac_up$histone <- ifelse(atac_up$gene %in% intersect(x$peaks,y$peaks), "H3K4me3+H3K27ac", atac_up$histone)
#pie3D(table(atac_up$histone),col = hcl.colors(4, "Earth"),labels=paste0(round(table(atac_up$histone)/length(atac_up$gene) * 100, 1), "%", ",", names(table(atac_up$histone))),explode = 0.05, theta=0.5, labelcex = 2)
pdf("~/Desktop/atac_2023_paper/DAR_up_histone_barplot.pdf")
ggplot(data=data.frame(table(atac_up$histone)), aes(x=Var1, y=Freq, fill=Var1)) +  geom_bar(stat="identity",color="black") + scale_fill_manual(values=c("#a6611a", "#dfc27d","#80cdc1","#018571")) + theme_minimal() + theme(text = element_text(size = 20),axis.text = element_text(size = 30),axis.text.x = element_blank(), axis.title.x = element_blank(), legend.text=element_blank(), legend.title=element_blank()) 
dev.off()
#overlap with DARs closed
atac_down <- atac_all[atac_all$logFC_treatment < 0,]
subsetByOverlaps(Signac::StringToGRanges(c(atac_down$gene)), unique(c(enhancer))) -> x
x$peaks <- paste(seqnames(x),ranges(x),sep="-")
subsetByOverlaps(Signac::StringToGRanges(c(atac_down$gene)), unique(c(promoter))) -> y
y$peaks <- paste(seqnames(y),ranges(y),sep="-")
atac_down$histone <- ifelse(atac_down$gene %in% x$peaks, "H3K27ac", ifelse(atac_down$gene %in% y$peaks, "H3K4me3", ifelse(atac_down$gene %in% intersect(x$peaks, y$peaks), "H3K4me3+H3K27ac+", "Other")))
atac_down$histone <- ifelse(atac_down$gene %in% intersect(x$peaks,y$peaks), "H3K4me3+H3K27ac+", atac_down$histone)
#pie3D(table(atac_down$histone),col = hcl.colors(4, "Earth"),labels=paste0(round(table(atac_down$histone)/length(atac_down$gene) * 100, 1), "%", ",", names(table(atac_down$histone))),explode = 0.05, theta=0.5, labelcex = 2)
pdf("~/Desktop/atac_2023_paper/DAR_down_histone_barplot.pdf")
ggplot(data=data.frame(table(atac_down$histone)), aes(x=Var1, y=Freq, fill=Var1)) +  geom_bar(stat="identity",color="black") + scale_fill_manual(values=c("#a6611a", "#dfc27d","#80cdc1","#018571")) + theme_minimal() + theme(text = element_text(size = 20),axis.text = element_text(size = 30),axis.text.x = element_blank(), axis.title.x = element_blank(), legend.text=element_blank(), legend.title=element_blank()) 
dev.off()

#overlap with cluster marker peaks
total_df <- data.frame()
path="~/List_of_SubClusters_Marker_Peaks/"
file_list <- list.files(path=path)
for(file_number in 1:length(file_list)){
  print(paste0(file_list[file_number]))
  #print(i)
  tmp_cluster <- read.csv(paste0(path,file_list[file_number]), sep="\t",header=TRUE)
  tmp_cluster$gene <- paste(tmp_cluster$seqnames, tmp_cluster$start, tmp_cluster$end, sep="-")
  tmp_cluster$celltype <- do.call(rbind,stringr::str_split(file_list[file_number],"_"))[,1]
  total_df <- rbind(total_df, tmp_cluster)
  }

subsetByOverlaps(Signac::StringToGRanges(total_df$gene), unique(enhancer)) -> x
x$peaks <- paste(seqnames(x),ranges(x),sep="-")
subsetByOverlaps(Signac::StringToGRanges(total_df$gene), unique(promoter)) -> y
y$peaks <- paste(seqnames(y),ranges(y),sep="-")
total_df$histone <- ifelse(total_df$gene %in% x$peaks, "H3K27ac", ifelse(total_df$gene %in% y$peaks, "H3K4me3", ifelse(total_df$gene %in% intersect(x$peaks, y$peaks), "H3K27ac+H3K4me3", "other")))
total_df$histone <- ifelse(total_df$gene %in% intersect(x$peaks,y$peaks), "H3K27ac+H3K4me3", total_df$histone)
pdf("~/scratch/cluster_marker_peak_histone_overlap.pdf",width=10)
pie3D(table(total_df$histone),col = hcl.colors(4, "Earth"),labels=paste0(round(table(total_df$histone)/length(total_df$gene) * 100, 1), "%", ",", names(table(total_df$histone))),explode = 0.05, theta=0.5, labelcex = 2)
dev.off()

#overlap with peaks in significant peak-to-gene linkages 
sub_p2g <- readRDS("~/Desktop/Subtype_p2g_redone.rds")
sub_p2g  <- sub_p2g %>% data.frame() %>% drop_na() %>% filter(Correlation > 0.45 & FDR < 1e-4)
#length(subsetByOverlaps(Signac::StringToGRanges(unique(sub_p2g$peakName)), unique(c(enhancer)))) / length(Signac::StringToGRanges(unique(sub_p2g$peakName)))
#length(subsetByOverlaps(Signac::StringToGRanges(unique(sub_p2g$peakName)), unique(c(promoter)))) / length(Signac::StringToGRanges(unique(sub_p2g$peakName)))

pdf("~/sub_p2g_histone.pdf",width=10)
subsetByOverlaps(Signac::StringToGRanges(unique(sub_p2g$peakName)), unique(c(enhancer))) -> x
x$peaks <- paste(seqnames(x),ranges(x),sep="-")
subsetByOverlaps(Signac::StringToGRanges(unique(sub_p2g$peakName)), unique(c(promoter))) -> y
y$peaks <- paste(seqnames(y),ranges(y),sep="-")
sub_p2g$histone <- ifelse(sub_p2g$peakName %in% x$peaks, "H3K27ac", ifelse(sub_p2g$gene %in% y$peaks, "H3K4me3", ifelse(sub_p2g$gene %in% intersect(x$peaks, y$peaks), "H3K4me3+H3K27ac+", "Other")))
sub_p2g$histone <- ifelse(sub_p2g$peakName %in% intersect(x$peaks,y$peaks), "H3K4me3+H3K27ac+", sub_p2g$histone)
pie3D(table(sub_p2g$histone),col = hcl.colors(4, "Earth"),labels=paste0(round(table(sub_p2g$histone)/length(sub_p2g$gene) * 100, 1), "%", ",", names(table(sub_p2g$histone))),explode = 0.05, theta=0.5, labelcex = 2)
dev.off()
#####################################################################################




