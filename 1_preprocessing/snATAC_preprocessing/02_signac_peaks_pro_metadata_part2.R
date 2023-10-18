library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(Matrix)

####pass argumets to R##
args <- commandArgs(TRUE)
subject <- noquote(args[1])
print(args)

Dir="/home/anjali5/scratch/graham_atac_analysis/"

#####cells from atac-pro raw_matrix######
cells <- read.csv(paste0(Dir,subject,"/raw_matrix/MACS2/barcodes.txt"), header=FALSE)
head(cells)
dim(cells)
####load peaks and fragment obj object#########
frags <- CreateFragmentObject(
path = paste0(Dir,subject,"/summary/fragments.tsv.gz"),
cells = cells$V1
)
print(frags)
peaks <- read.table(paste0(Dir,subject,"/peaks/MACS2/peaks.bed"),sep="\t",header=FALSE, col.names = c("chr", "start", "end", "width","strand","name","score", "fold_change","neg_log10pvalue_summit","neg_log10qvalue_summit","relative_summit_position"))
head(peaks)
dim(peaks)
#####convert peaks to GRanges#####
peaks <- makeGRangesFromDataFrame(peaks)
head(peaks)
length(peaks)
####metadata from atac-pro output####
md <- read.table(
file = paste0(Dir,subject,"/summary/atac.MACS2.qc_per_barcode.txt"),
stringsAsFactors = FALSE,
sep = "\t",
header = TRUE
)
dim(md)
#####cells from atac-pro raw_matrix######
#cells <- read.csv(paste0(Dir,subject,"/raw_matrix/MACS2/barcodes.txt"), header=FALSE)
#head(cells)
######filter metadata from atac-pro######
md.filter <- subset(md, md$bc %in% cells$V1)
rownames(md.filter) <- md.filter[,1]
dim(md.filter)
####metadata from 10x output####
library(data.table)
md10x <- read.table(
file = paste0(Dir,subject,"/singlecell.csv"),
stringsAsFactors = FALSE,
row.names = 1, sep = ",",
header = TRUE
)[-1, ] 
#md10x$org_cell <- md10x$is__cell_barcode
#md10x$is__cell_barcode <- ifelse(md10x$is__cell_barcode=="0", "1", "1") 
#write.table(md10x,paste0(Dir, "atac_", subject, "/singlecell.csv"), quote=F, row.names=F, col.names = F, sep = ',')
#md10x$barcode <- iflese(md$excluded_reason=="1" | md$excluded_reason=="3", paste0(md$barcode,"-E"),md$bar>
print(head(md10x))
print(dim(md10x))
extra <- dplyr::filter(md.filter,!(md.filter$bc %in% rownames(md10x)))
dim(extra)
print("no extra if 0 rows  above")
#####filter 10x metadata######
md10x.filter <- subset(md10x, rownames(md10x) %in% cells$V1)
md10x.filter$bc <- rownames(md10x.filter)
print(dim(md10x.filter))
#####make sure md and md10x match#####
nrow(md.filter) == nrow(md10x.filter)
metadata = merge(md.filter, md10x.filter, by="bc")
rownames(metadata) <- metadata$bc
head(metadata)
dim(metadata)
#######make a seurat object#########
counts <- FeatureMatrix(
fragments = frags,
features = peaks,
cells = cells$V1
)
print(head(counts))
print(dim(counts))
####chromatin_assay###
assay <- CreateChromatinAssay(counts, fragments = frags)
print(assay)
#####Seurat_obj####
atac <- CreateSeuratObject(
counts = assay,
assay = 'peaks',
project = 'ATAC',
meta.data = metadata
)
print(atac)
print(head(atac@meta.data))
##SAVE 
saveRDS(atac,paste0(Dir,subject,"/atac.rds"))
print("obj done and saved")

