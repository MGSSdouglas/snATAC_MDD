library(ArchR)
library(readr)
library(stringr)
library(GenomicRanges)
library(Signac)


#Load in all input data

input_dir <- "~/projects/rrg-gturecki/Anjali_snATAC/ArchR/"
p2g_input_dir <- "~/projects/rrg-gturecki/Anjali_snATAC/snps_to_plot/peak_to_gene_linkages/"

#proj <- loadArchRProject(
    #paste0(input_dir, "ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered"))

#broad_proj <- loadArchRProject(
    #paste0(input_dir, "ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered"))


#Load in all the information needed for visualization:

source("getHiC.r")
source("plot_atac_connections.r")

#Note - to get all the most recent ArchR functions required for this code to work, clone
#the most recent version of ArchR from Github and load in all R scripts
#In the present working directory - 
#git clone "https://github.com/GreenleafLab/ArchR/"

for(file in list.files("./ArchR/R/")) {
     source(paste0("./ArchR/R/",file))
}
source("ArchRBrowserWithExternalData.r")

#corCutOff <- 0.45
corCutOff <- 0.45
#Set up the following GRanges objects for future use:
#broad_p2g, fine_p2g, broad_p2p, fine_p2p

#Get fine p2p links, including negative correlations
fine_peakToPeakDF <- metadata(proj@peakSet)$CoAccessibility
#fine_peakToPeakDF <- fine_peakToPeakDF[which(abs(fine_peakToPeakDF$correlation) > corCutOff),]
fine_peakToPeakDF <- fine_peakToPeakDF[which(fine_peakToPeakDF$correlation > corCutOff),]
fine_peakToPeakDF$peak1 <- (metadata(fine_peakToPeakDF)$peakSet %>% {paste0(seqnames(.), "-", start(.), "-", end(.))})[fine_peakToPeakDF$queryHits]
fine_peakToPeakDF$peak2 <- (metadata(fine_peakToPeakDF)$peakSet %>% {paste0(seqnames(.), "-", start(.), "-", end(.))})[fine_peakToPeakDF$subjectHits]

fine_p2p_peak1_granges <- (metadata(fine_peakToPeakDF)$peakSet)[fine_peakToPeakDF$queryHits]
fine_p2p_peak2_granges <- (metadata(fine_peakToPeakDF)$peakSet)[fine_peakToPeakDF$subjectHits]
fine_peak1_midpoints <- 0.5*(start(fine_p2p_peak1_granges)+end(fine_p2p_peak1_granges))
fine_peak2_midpoints <- 0.5*(start(fine_p2p_peak2_granges)+end(fine_p2p_peak2_granges))
fine_p2p_coordinates <- cbind(fine_peak1_midpoints, fine_peak2_midpoints)
fine_p2p_coordinates <-  t(apply(fine_p2p_coordinates, 1, sort)) #this will take a few minutes

fine_p2p <- GRanges(data.frame(chr=seqnames(fine_p2p_peak1_granges),
                               start=fine_p2p_coordinates[,1],
                               end=fine_p2p_coordinates[,2],
                               correlation=fine_peakToPeakDF$correlation,
                               FDR=fine_peakToPeakDF$FDR,
                               peak1=fine_peakToPeakDF$peak1,
                               peak2=fine_peakToPeakDF$peak2))
#if(FALSE){
#Repeat for broad p2p
broad_peakToPeakDF <- metadata(broad_proj@peakSet)$CoAccessibility
#broad_peakToPeakDF <- broad_peakToPeakDF[which(abs(broad_peakToPeakDF$correlation) > corCutOff),]
broad_peakToPeakDF <- broad_peakToPeakDF[which(broad_peakToPeakDF$correlation) > corCutOff,]
broad_peakToPeakDF$peak1 <- (metadata(broad_peakToPeakDF)$peakSet %>% {paste0(seqnames(.), "-", start(.), "-", end(.))})[broad_peakToPeakDF$queryHits]
broad_peakToPeakDF$peak2 <- (metadata(broad_peakToPeakDF)$peakSet %>% {paste0(seqnames(.), "-", start(.), "-", end(.))})[broad_peakToPeakDF$subjectHits]

broad_p2p_peak1_granges <- (metadata(broad_peakToPeakDF)$peakSet)[broad_peakToPeakDF$queryHits]
broad_p2p_peak2_granges <- (metadata(broad_peakToPeakDF)$peakSet)[broad_peakToPeakDF$subjectHits]
broad_peak1_midpoints <- 0.5*(start(broad_p2p_peak1_granges)+end(broad_p2p_peak1_granges))
broad_peak2_midpoints <- 0.5*(start(broad_p2p_peak2_granges)+end(broad_p2p_peak2_granges))
broad_p2p_coordinates <- cbind(broad_peak1_midpoints, broad_peak2_midpoints)
broad_p2p_coordinates <-  t(apply(broad_p2p_coordinates, 1, sort))

broad_p2p <- GRanges(data.frame(chr=seqnames(broad_p2p_peak1_granges),
                               start=broad_p2p_coordinates[,1],
                               end=broad_p2p_coordinates[,2],
                               correlation=broad_peakToPeakDF$correlation,
                               FDR=broad_peakToPeakDF$FDR,
                               peak1=broad_peakToPeakDF$peak1,
                               peak2=broad_peakToPeakDF$peak2))
#}
#Get peak to gene coordinates, including negative correlations and gene TSS,
#Get all gene TSS values
#fine_tss_coords <- getTSS(proj) IGNORE DOING THIS WITH TSS FOR NOW. The output of getTSS isn't the same length
#as  metadata(fine_p2geneDF)$geneSet), and it isn't obvious how to get a dictionary going from gene --> TSS

#we will use start(metadata(fine_p2geneDF)$geneSet) since every gene corresponds to a single position (not a range)
#looking this up for a few examples, it isn't obvious where the position falls within the gene, and how this is computed
#it is within the gene but not necessarily the TSS - we can figure this out later

fine_p2geneDF <- readRDS(paste0(p2g_input_dir, "Subtype_p2g_redone.rds"))
#fine_p2geneDF <- fine_p2geneDF[which(abs(fine_p2geneDF$Correlation) > corCutOff),]
fine_p2geneDF <- fine_p2geneDF[which(fine_p2geneDF$Correlation) > corCutOff,]

fine_p2geneDF$peakName <- (metadata(fine_p2geneDF)$peakSet %>% {paste0(seqnames(.), "-", start(.), "-", end(.))})[fine_p2geneDF$idxATAC]
fine_p2geneDF$geneName <- mcols(metadata(fine_p2geneDF)$geneSet)$name[fine_p2geneDF$idxRNA]
#fine_p2geneDF$TSS <- start(fine_tss_coords)[fine_p2geneDF$idxRNA]
fine_p2geneDF$gene_connection_pos <- start(metadata(fine_p2geneDF)$geneSet)[fine_p2geneDF$idxRNA]
fine_p2g_peak <- (metadata(fine_p2geneDF)$peakSet)[fine_p2geneDF$idxATAC]
fine_p2g_peak_midpoints <- 0.5*(start(fine_p2g_peak)+end(fine_p2g_peak))
fine_p2g_coordinates <- cbind(fine_p2g_peak_midpoints, fine_p2geneDF$gene_connection_pos)
fine_p2g_coordinates <-  t(apply(fine_p2g_coordinates, 1, sort))
fine_p2g <- GRanges(data.frame(chr=seqnames(fine_p2g_peak),
                               start=fine_p2g_coordinates[,1],
                               end=fine_p2g_coordinates[,2],
                               correlation=fine_p2geneDF$Correlation,
                               FDR=fine_p2geneDF$FDR,
                               peak=fine_p2geneDF$peakName,
                               peakType=fine_p2geneDF$peakType,
                               gene=fine_p2geneDF$geneName))
                               #TSS=fine_p2geneDF$TSS ))
#can always add more info here as needed
#for comparison:
#fine_p2g_from_proj <- getPeak2GeneLinks(proj, corCutOff=corCutOff)$Peak2GeneLinks
#if(FALSE){
#Repeat for broad project
#broad_tss_coords <- getTSS(broad_proj)
broad_p2geneDF <- readRDS(paste0(p2g_input_dir, "broad_p2g.rds"))
#broad_p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
#broad_p2geneDF <- broad_p2geneDF[which(abs(broad_p2geneDF$Correlation) > corCutOff),]
broad_p2geneDF <- broad_p2geneDF[which(broad_p2geneDF$Correlation > corCutOff),]

broad_p2geneDF$peakName <- (metadata(broad_p2geneDF)$peakSet %>% {paste0(seqnames(.), "-", start(.), "-", end(.))})[broad_p2geneDF$idxATAC]
broad_p2geneDF$geneName <- mcols(metadata(broad_p2geneDF)$geneSet)$name[broad_p2geneDF$idxRNA]
#broad_p2geneDF$TSS <- start(broad_tss_coords)[broad_p2geneDF$idxRNA]
broad_p2geneDF$gene_connection_pos <- start(metadata(broad_p2geneDF)$geneSet)[broad_p2geneDF$idxRNA]
broad_p2g_peak <- (metadata(broad_p2geneDF)$peakSet)[broad_p2geneDF$idxATAC]
broad_p2g_peak_midpoints <- 0.5*(start(broad_p2g_peak)+end(broad_p2g_peak))
broad_p2g_coordinates <- cbind(broad_p2g_peak_midpoints, broad_p2geneDF$gene_connection_pos)
broad_p2g_coordinates <-  t(apply(broad_p2g_coordinates, 1, sort))
broad_p2g <- GRanges(data.frame(chr=seqnames(broad_p2g_peak),
                               start=broad_p2g_coordinates[,1],
                               end=broad_p2g_coordinates[,2],
                               correlation=broad_p2geneDF$Correlation,
                               FDR=broad_p2geneDF$FDR,
                               peak=broad_p2geneDF$peakName,
                               peakType=broad_p2geneDF$peakType,
                               gene=broad_p2geneDF$geneName))
                               #TSS=broad_p2geneDF$TSS ))
#}
#Write these four GRanges objects as outputs to be reused

output_dir <- "~/projects/rrg-gturecki/snp_visualization/data/"
write.table(fine_p2p, paste0(output_dir,"fine_p2p_", corCutOff, ".tsv"),
                          row.names=FALSE, sep="\t", quote=FALSE)

write.table(fine_p2g, paste0(output_dir,"fine_p2g_", corCutOff, ".tsv"),
                          row.names=FALSE, sep="\t", quote=FALSE)

write.table(broad_p2p, paste0(output_dir,"broad_p2p_", corCutOff, ".tsv"),
                          row.names=FALSE, sep="\t", quote=FALSE)

write.table(broad_p2g, paste0(output_dir,"broad_p2g_", corCutOff, ".tsv"),
                          row.names=FALSE, sep="\t", quote=FALSE)


