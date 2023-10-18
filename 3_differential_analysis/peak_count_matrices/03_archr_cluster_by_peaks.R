
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Matrix)
set.seed(1234)
library(future)
library(tidyverse)
#library(ArchR)
cpu_n = 40

options(future.globals.maxSize = 220000 * 1024 ^ 2)
plan("multiprocess", workers = cpu_n)
plan()

#load concatenated count data
counts <- readMM("/home/anjali5/projects/def-cnagy/For_Wenmin/For_Doruk/merged_case_control.v4/all_adata_counts.005.mtx")
counts <- t(counts)
head(counts)
dim(counts)
barcodes <- read_delim("/home/anjali5/projects/def-cnagy/For_Wenmin/For_Doruk/merged_case_control.v4/all_adata_obs.005.tsv", delim="\t")
head(barcodes)
dim(barcodes)
peaks <- read_delim("/home/anjali5/projects/def-cnagy/For_Wenmin/For_Doruk/merged_case_control.v4/all_adata_var.005.tsv", delim="\t")
head(peaks)
dim(peaks)

colnames(counts) <- barcodes$x
rownames(counts) <- peaks$x
rownames(barcodes) <- barcodes$x

head(counts)
dim(counts)

#create peak assay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = NULL,
  fragments = NULL,
  min.cells = 0,
  min.features = 0
)
chrom_assay
brain  <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = data.frame(barcodes)
)
brain
#noramlize signac object
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(object = brain)

#cluster signac object
brain <- RunUMAP(object = brain, reduction = 'lsi', dims = 2:40)
brain <- FindNeighbors(object = brain, reduction = 'lsi', dims = 2:40)
brain <- FindClusters(object = brain, verbose = FALSE, algorithm = 3,resolution = 2.5)

brain

saveRDS(brain, "~/scratch/ArchR_cluster_by_peaks.rds")

#add genome information 
chrom_hg38 <- read.table("~/scratch/graham_merge_analysis/R_scripts/chromInfo.txt")
annotations <- readRDS("/home/anjali5/scratch/UCSC_annotations.rds")
#annotations
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- chrom_hg38

# add the gene information to the object
Annotation(brain) <- annotations
gene.activities <- GeneActivity(brain)

# add the gene activity matrix to the Seurat object as a new assay
brain[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)

saveRDS(brain, "~/scratch/ArchR_cluster_by_peaks_updated.rds")

