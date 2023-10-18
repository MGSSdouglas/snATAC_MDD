library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(Matrix)

Dir="/home/anjali5/scratch/graham_atac_analysis/"
####pass argumets to R##
args <- commandArgs(TRUE)
subject <- noquote(args[1])
print(args)
####load seurat object#########
atac <- readRDS(paste0(Dir,"atac_", subject,"/atac.rds"))
atac
####load ArchR barcodes######
barcodes <- read.csv(paste0(Dir,"atac_", subject,"/archR_barcodes"),header=TRUE,sep="")
head(barcodes)
dim(barcodes)
#####amulet output##################
MultipletBarcodes <- read.csv(paste0(Dir, "atac_", subject, "/amulet_archR/MultipletBarcodes_01.txt"),header=FALSE)
head(MultipletBarcodes)
dim(MultipletBarcodes)
head(md)
dim(md)
#######Add Amulet doublets to metadata###
md <- atac@meta.data
md$amulet <- ifelse(md$bc %in% MultipletBarcodes$V1, "1", "0")
md$barcode <- md$bc
head(md)
#md.filter$bc <- ifelse(md.$bc %in% MultipletBarcodes$V1, paste0(md.filter$bc,"-D"), md.filter$bc)
########Save ArchR atac.object#############
atac <- AddMetaData(atac, md)
atac
saveRDS(atac,paste0(Dir,"atac_",subject,"/atac_updated.rds"))
atac_archR <- subset(atac, bc %in% barcodes$barcode)
atac_archR
saveRDS(atac_archR,paste0(Dir,"atac_",subject,"/atac_archR.rds"))
print("done")
###########################################################

