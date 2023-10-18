library(ArchR)
library(readr)
library(stringr)
library(GenomicRanges)
library(Signac)

#*****************************************************************************
#Load in all input data

input_dir <- "~/projects/rrg-gturecki/Anjali_snATAC/ArchR/"
proj <- loadArchRProject(
    paste0(input_dir, "ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered"))

broad_proj <- loadArchRProject(
    paste0(input_dir, "ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered"))


#Load in all the information needed for visualization:

source("getHiC.r")
source("plot_atac_connections.r")

#Note - to get all the most recent ArchR functions required for this code to work, clone
#the most recent version of ArchR from Github and load in all R scripts
#In the present working directory - 
#git clone "https://github.com/GreenleafLab/ArchR/"

for(file in list.files("./ArchR/R/")) { #assume working directory is vikram_visualizations
     source(paste0("./ArchR/R/",file))
}
source("ArchRBrowserWithExternalData.r")

corCutOffToLoad <- 0.45

#Load in p2p and p2g data
#OPTION - rather than loading in external peakToPeak/peakToGene, if we only want positive connections:

#broad_p2p_list <- getCoAccessibility(broad_proj, corCutOff=corCutOffToLoad)
#broad_p2g_list <- getPeak2GeneLinks(broad_proj, corCutOff=corCutOffToLoad)

#comment everything below out.
input_dir <- "~/projects/rrg-gturecki/vikram_visualization/data/"
broad_p2p <- GRanges(read_tsv(paste0(input_dir, "broad_p2p_", corCutOffToLoad, ".tsv")))
broad_p2p_list <- list()

fine_p2p <- GRanges(read_tsv(paste0(input_dir, "fine_p2p_", corCutOffToLoad, ".tsv")))
fine_p2p_list <- list()

input_df <- read_tsv(paste0(input_dir, "broad_p2g_", corCutOffToLoad, ".tsv"))
broad_p2g <- GRanges(input_df)
#broad_p2g <- GRanges(read_tsv(paste0(input_dir, "broad_p2g_", corCutOffToLoad, ".tsv")))

fine_p2g <- GRanges(read_tsv(paste0(input_dir, "fine_p2g_", corCutOffToLoad, ".tsv")))
fine_p2g_list <- list()

broad_p2g_list <- list()

#Renaming
colnames(mcols(fine_p2p))[1] <- "value"
colnames(mcols(fine_p2g))[1] <- "value"
colnames(mcols(broad_p2p))[1] <- "value"
colnames(mcols(broad_p2g))[1] <- "value"

#List creation
broad_p2p_list[["Broad P2P Co-Accessibility"]] <- broad_p2p
fine_p2p_list[["Fine P2P Co-Accessibility"]] <- fine_p2p
fine_p2g_list[["Fine P2G Linkage"]] <- fine_p2g
broad_p2g_list[["Broad P2G Linkage"]] <- broad_p2g


#Load in histone data
#corresponds to NeuN+/NeuN- and H3K27ac and H3K4me3
histone_dir <- "~/projects/rrg-gturecki/vikram_visualization/data/"
pos_me_df <- read.table(paste0(histone_dir, "NeuN+_H3K4me3_hg38.bed"))
pos_me_granges <- GRanges(data.frame(seqnames=pos_me_df[,1], start=pos_me_df[,2], end=pos_me_df[,3]))
names(pos_me_granges) <- rep("NeuN+ H3K4me3", length(pos_me_granges))

pos_ac_df <- read.table(paste0(histone_dir,"NeuN+_H3K27ac_hg38.bed"))
pos_ac_granges <- GRanges(data.frame(seqnames=pos_ac_df[,1], start=pos_ac_df[,2], end=pos_ac_df[,3]))
names(pos_ac_granges) <- rep("NeuN+ H3K27ac", length(pos_ac_granges))


neg_me_df <- read.table(paste0(histone_dir, "NeuN-_H3K4me3_hg38.bed"))
neg_me_granges <- GRanges(data.frame(seqnames=neg_me_df[,1], start=neg_me_df[,2], end=neg_me_df[,3]))
names(neg_me_granges) <- rep("NeuN- H3K4me3", length(neg_me_granges))

neg_ac_df <- read.table(paste0(histone_dir, "NeuN-_H3K27ac_hg38.bed"))
neg_ac_granges <- GRanges(data.frame(seqnames=neg_ac_df[,1], start=neg_ac_df[,2], end=neg_ac_df[,3]))
names(neg_ac_granges) <- rep("NeuN- H3K27ac", length(neg_ac_granges))

#Convert histone data into a GRangesList for visualization
histone_list <- GRangesList(pos_me_granges, pos_ac_granges, neg_me_granges, neg_ac_granges)
names(histone_list) <- c("NeuN+ H3K4me3","NeuN+ H3K27ac", "NeuN- H3K4me3", "NeuN- H3K27ac")
region_to_assess <- GRanges(data.frame(seqnames="chr2", start=234440000, end=234500000))
for(i in c(1:4)) {
    print(subsetByOverlaps(histone_list[[i]], region_to_assess)) #is there only one type of histone modification present here??
}

#Enter in SNPs for visualization - make each one a separate peak track

#Load in HiC fragments/loops
hiC_dir <- "~/projects/rrg-gturecki/vikram_visualization/data/"
neu_neg_path <-  paste0(hiC_dir,"neu_neg_loops_hg38.bed")
neu_pos_path <-  paste0(hiC_dir, "neu_pos_loops_hg38.bed")

all_hiC_loops <- get_hiC_loops(neu_neg_path, neu_pos_path)
all_hiC_fragments <- get_hiC_fragments(neu_neg_path, neu_pos_path)

#*******************************************************************************
#Visualize.

#row 2 of excel sheet
snp <- GRanges(data.frame(seqnames="chr2", start="234446420", end="234446421"))
LDsnp <- GRanges(data.frame(seqnames="chr2", start="23444600", end="234446601"))

snp_list <- GRangesList(snp, LDsnp)
names(snp_list) <- c("sSNP", "LD SNPs")
snp_extend <- 100
snp_extended <- snp
start(snp_extended) <- start(snp_extended) - snp_extend
end(snp_extended) <- end(snp_extended) + snp_extend

peak <- GRanges(data.frame(seqnames="chr2", start="234446364", end="234446864"))
arl4_tss <- GRanges(data.frame(seqnames="chr2", start="234493041", end="234493042"))

external_features <- list(snp_extended, histone_list, all_hiC_fragments)
external_feature_names <- c("sSNP", "Histone ChIP-seq", "HiC Fragments")
p2p_list <- broad_p2p_list
external_loops <- c(p2p_list, all_hiC_loops)
external_loop_names <- c("Peak-To-Peak Co-Accessibility", "HiC Connections")

broad_proj$subcluster_condition <- paste(broad_proj$SubClusters, broad_proj$condition, sep="_")
proj$subcluster_condition <- paste(proj$SubClusters, proj$condition, sep="_")

plot_ATAC_connections_two_highlighted_regions(proj = broad_proj, left_extension=10000, right_extension=10000, highlight_extend=1000,
                    region_1=peak, region_2=arl4_tss, chr="chr2", groupBy="subcluster_condition",
                    useGroups=c("Mic2_Case", "Mic2_Control"),#new!!
                    baseSize=6, facetbaseSize=6,
                    external_features=external_features, external_feature_names=external_feature_names,
                    peak_to_gene_loops=broad_p2g_list, 
                    external_loops=external_loops, external_loop_names=external_loop_names,
                    prefix="rs10210757_subcluster")


#Row 3 of excel sheet - rs5995992
snp <- GRanges(data.frame(seqnames="chr22", start=41091214, end=41091215))
snp_extend <- 100
snp_extended <- snp
start(snp_extended) <- start(snp_extended) - snp_extend
end(snp_extended) <- end(snp_extended) + snp_extend

peak <- GRanges(data.frame(seqnames="chr22", start="41091360", end="41091860"))
CSDC2_tss <- GRanges(data.frame(seqnames="chr22", start="41561010", end="41561011"))
TOB2_tss <- GRanges(data.frame(seqnames="chr22", start="41433494", end="41433494"))
external_features <- list(snp_extended, histone_list, all_hiC_fragments)
proj$subcluster_sex_condition <- paste(proj$SubClusters, proj$sex, proj$condition, sep="_")
external_feature_names <- c("sSNP", "Histone", "HiC Fragments")
p2p_list <- fine_p2p_list
external_loops <- c(p2p_list, all_hiC_loops)
external_loop_names <- c("Peak-To-Peak Co-Accessibility", "HiC Connections")

plot_ATAC_connections_two_highlighted_regions(proj = proj, left_extension=10000, right_extension=10000, highlight_extend=1000,
                    peak, TOB2_tss, chr="chr22", groupBy="subcluster_sex_condition",
                    useGroups=c("Mic2_male_Case", "Mic2_male_Control", "Mic2_female_Case", "Mic2_female_Control" ),#new!!
                    baseSize=6, facetbaseSize=6,
                    external_features=external_features, external_feature_names=external_feature_names,
                    peak_to_gene_loops=fine_p2g_list, 
                    external_loops=external_loops, external_loop_names=external_loop_names,
                    prefix="rs5995992_TOB2")

#Row 7 of excel sheet: rs80079217
snp <- GRanges(data.frame(seqnames="chr2", start=160363479, end=160363480))
snp_extend <- 100
snp_extended <- snp
start(snp_extended) <- start(snp_extended) - snp_extend
end(snp_extended) <- end(snp_extended) + snp_extend

peak <- GRanges(data.frame(seqnames="chr2", start="160363031", end="160363531"))
MIR4785_tss <- GRanges(data.frame(seqnames="chr2", start="160407810", end="160407811"))
external_features <- list(snp_extended, histone_list, all_hiC_fragments)
proj$subcluster_condition <- paste(proj$SubClusters, proj$condition, sep="_")
external_feature_names <- c("sSNP", "Histone", "HiC Fragments")
p2p_list <- fine_p2p_list
external_loops <- c(p2p_list, all_hiC_loops)
external_loop_names <- c("Peak-To-Peak Co-Accessibility", "HiC Connections")

plot_ATAC_connections_two_highlighted_regions(proj = proj, left_extension=10000, right_extension=10000, highlight_extend=1000,
                    peak, MIR4785_tss, chr="chr2", groupBy="subcluster_condition",
                    useGroups=c("ExN1_Case", "ExN1_Control"),#new!!
                    baseSize=6, facetbaseSize=6,
                    external_features=external_features, external_feature_names=external_feature_names,
                    peak_to_gene_loops=fine_p2g_list, 
                    external_loops=external_loops, external_loop_names=external_loop_names,
                    prefix="rs80079217_MIR4785")

#Row 8 of excel sheet: rs174561 
snp <- GRanges(data.frame(seqnames="chr11", start=61815236, end=61815237))
snp_extend <- 100
snp_extended <- snp
start(snp_extended) <- start(snp_extended) - snp_extend
end(snp_extended) <- end(snp_extended) + snp_extend

peak <- GRanges(data.frame(seqnames="chr11", start="61815147", end="61815647"))
FADS1_tss <- GRanges(data.frame(seqnames="chr11", start="61799627", end="61799628"))
TMEM258_tss <- GRanges(data.frame(seqnames="chr11", start="61768501", end="61768502"))
external_features <- list(snp_extended, histone_list, all_hiC_fragments)
proj$subcluster_condition <- paste(proj$SubClusters, proj$condition, sep="_")
external_feature_names <- c("sSNP", "Histone", "HiC Fragments")
p2p_list <- fine_p2p_list
external_loops <- c(p2p_list, all_hiC_loops)
external_loop_names <- c("Peak-To-Peak Co-Accessibility", "HiC Connections")

plot_ATAC_connections_two_highlighted_regions(proj = proj, left_extension=10000, right_extension=10000, highlight_extend=1000,
                    peak, FADS1_tss, chr="chr11", groupBy="SubClusters",
                    useGroups=c("Ast1", "Ast2", "Ast3", "Ast4"),#new!!
                    baseSize=6, facetbaseSize=6,
                    external_features=external_features, external_feature_names=external_feature_names,
                    peak_to_gene_loops=fine_p2g_list, 
                    external_loops=external_loops, external_loop_names=external_loop_names,
                    prefix="rs174561_FADS1")


#Default: everything.
#bulk track - make the size more reasonable
#gene track
#HiC - same coloring as
#Histone - colored by neu+/neu-/me/ac
#P2P - only connections to SNP (like before)
#P2G - leave everything
#sSNP

#Get rid of default "LoopTrack" option when nothing is found in the window

#Create groups for cases/controls (ex: Mic1_case, Mic1_control within same project use that as the groupBy)
#groupBy="group"
#compare getPeak2GeneLinks output with our custom p2g links with TSS that also include negative



