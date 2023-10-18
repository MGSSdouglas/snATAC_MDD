library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(Matrix)
library(muscat)
library(DESeq2)
#library(pcaExplorer)
library(tidyverse)
#library(PCAtools)
library(limma)
library(scater)
library(BiocParallel)
library(edgeR)
library(future)
library(ArchR)
library(dplyr)

cpu_n = 40

options(future.globals.maxSize = 220000 * 1024 ^ 2)
plan("multiprocess", workers = cpu_n)
plan()

#load broad and cluster peak matrices in seurat object
obj <- readRDS("/home/anjali5/projects/def-cnagy/anjali5/seuratObj_Extended_Broad_PeakMatrix_ArchR.rds") #broad

#subset subject with misassigned reported sex
obj  <- subset(obj, ident=c("atac_A32B32female"),invert=TRUE)
obj
print(table(obj$Subject,obj$sex))
print(table(obj$Subject,obj$condition))

#convert Seurat object to single-cell experiment class
Idents(obj) <- "ClustersMapped" 
sce <- as.SingleCellExperiment(obj)
sce

#filter low-quality peaks
#sce <- sce[rowSums(counts(sce) > 0) > 0, ]
#saveRDS(sce, "~/projects/def-gturecki/anjali5/Finalized_DAR/subcluster/sce.rds")

#Pseudobulk per subject per cluster using muscat
sce <- prepSCE(sce,kid = "ClustersMapped",gid = "condition",sid = "Subject",drop=FALSE)
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids
t(table(sce$cluster_id, sce$sample_id))
#saveRDS(sce,"~/scratch/sce_broad.rds")

##aggregate peak counts per subjects per cluster##########################
pb <- aggregateData(sce,assay = "counts", fun = "sum",by = c("cluster_id","sample_id"))
assayNames(pb)
t(head(assay(pb)))
#pb
#filter out chrY peaks
pb <- pb[-(grep("^chrY:", row.names(pb))),]
saveRDS(pb, "~/projects/def-gturecki/anjali5/Finalized_DAR/subcluster/pb_broad.rds")
#pb <- readRDS("~/projects/def-gturecki/anjali5/Finalized_DAR/subcluster/pb_broad.rds")
#pb
#check_samples
#pb <- pb[,colData(pb)$BrainID %ni% check_samples] #sensitivity check
#print(pb)

##get per-subject metadata
temp <- colData(pb)
colData(pb) <- temp
pb
temp$sample_id <- rownames(temp)
table(temp$sample_id,temp$group_id)

##match aggregated count metadata with single-cell metadata
#get cell counts per subject to add as covar
#get metdata from object 
metadata <- obj@meta.dat
NumCells <- metadata %>% group_by(Subject) %>% summarise(n = n())
NumCells <- NumCells[match(temp$sample_id, NumCells$Subject),]
temp$NumCells <- NumCells$n

##create design matrix
anno <- unique.data.frame(temp[c("sample_id", "group_id", "sex", "batch", "Age", "pHValue", "PMI", "NumCells")])
anno$PMI <- as.numeric(anno$PMI)
anno$Age <- as.numeric(anno$Age)
anno$pHValue <- as.numeric(anno$pHValue)
anno$batch <- as.factor(anno$batch)
anno$sex <- as.factor(anno$sex)
#head(anno)
anno$group2 <- paste(anno$group_id, anno$sex, sep = ".")
design <- model.matrix(~ 0 + group2 + batch + Age + PMI +NumCells, data = anno)
dimnames(design) <- list(anno$sample_id, levels = make.names(colnames(design)))
#head(design)

#make contrasts & run model
cm <- makeContrasts(treatment = (group2case.male + group2case.female) / 2 - (group2control.male+ group2control.female) / 2,levels = make.names(colnames(design)))
result <- pbDS(pb, design = design, contrast = cm, method = "limma-voom", min_cells=10, filter = "both") #min_cells=10 for broad cell-type 

#save results
saveRDS(result, "~/scratch/broad_result_limma_allPeaks_corrected_filtered.rds")

##Matt's method (gives equivalent result)
#anno$group <- ifelse(anno$group_id=="case",1,0)
#anno$sex <- ifelse(anno$sex=="male",1,0)
#design <- model.matrix(~ 0 + group + sex + sex*group + batch + Age + PMI + NumCells, data = anno)
#cm <- makeContrasts(treatment = group+ (group.sex)/2, 
#                    male.trt = group+group.sex, female.trt = group,
#                    interaction = group.sex,
#                    levels=make.names(colnames(design)))
#rm(obj)
#result2 <- pbDS(pb, design = design, contrast = cm, method = "limma-voom", min_cells=5, filter = "both")

##for each subject and group get fraction of cells with accessible peak frequencies 
frq <- calcExprFreqs(sce, assay = "counts", th = 0)
head(frq)
head(assay(frq))

# expression frequencies in each
# sample & group; 1st cluster
t(head(assay(frq), 5))
#saveRDS(frq, "~/scratch/frq_broad.rds")
#coefs <- c("treatment","male.trt", "female.trt", "interaction")
res_freq <- resDS.mod(sce, result, cpm = FALSE, frq = frq, bind = "col")
saveRDS(res_freq, "~/scratch/res_freq_broad.rds")

#filter results to keep peaks accessible in an average of 1% of cells in at least 1 group
#Load result file
atac <- readRDS("~/projects/def-gturecki/anjali5/Finalized_DAR/subcluster/res_freq_broad.rds")
tmp   <- atac[, grep("case.frq|control.frq", colnames(atac))]
tmp$Greater1pct <- ifelse(tmp$case.frq  > 0.01 | tmp$control.frq > 0.01, "TRUE", "FALSE")
atac$Greater1pct <- tmp$Greater1pct

#write output to file
atac_sub <- atac[abs(atac$logFC_treatment) > log2(1.1) & atac$p_adj.loc_treatment < 0.05 & atac$Greater1pct==TRUE,] 
table(atac_sub$cluster_id)
write.csv(atac_sub, "~/scratch/DARs_broad.csv")

print("Done")
sessionInfo()



