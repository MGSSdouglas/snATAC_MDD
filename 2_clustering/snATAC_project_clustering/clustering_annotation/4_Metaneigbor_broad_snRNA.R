#MetaNeighbor

library(MetaNeighbor)
library(Seurat)
library(SingleCellExperiment)
library(ArchR)
#Use Allen Brain motor cortex referece - using human M1 cortex data from https://portal.brain-map.org/atlases-and-data/rnaseq

library(tidyverse)
sessionInfo()
traceback()

print("Loading snATAC data")

ArchR <- readRDS("~/projects/def-gturecki/anjali5/ArchR_Subcluster_by_peaks_updated.rds")
ArchR$SubClusters <- plyr::mapvalues(ArchR$SubClusters, from = unique(ArchR$SubClusters), to = c("ExN12","ExN9", "Mix2","ExN2", "ExN6", "ExN10","ExN8", "ExN11", "ExN1","InPV", "Oli7", "Ast2", "Oli5","ExN4", "ExN5","ExN3", "InLAMP5","InSST", "OPC2", "InVIP", "Oli1", "Oli4", "End1", "Ast3", "OPC4", "End2", "Mic1", "Mic2", "Oli3", "Mix1", "Oli2", "Ast4", "ExN7", "OPC1", "OPC3", "InN3", "OLi6", "Ast1"))

DefaultAssay(ArchR) <- "RNA"

ArchR <- FindVariableFeatures(ArchR, selection.method = "vst", nfeatures = 3000)
ArchR@meta.data["study_id"] <- rep("snATAC", dim(ArchR)[2])
ArchR_var_genes <- VariableFeatures(ArchR)

#https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
ABI_metadata <- read_csv("/home/anjali5/projects/def-cnagy/anjali5/ABI/metadata.csv")
ABI_metadata <- column_to_rownames(ABI_metadata, "sample_name")

ABI_data <- read_csv("/home/anjali5/projects/def-cnagy/anjali5/ABI/matrix.csv")
ABI_data <- column_to_rownames(ABI_data, "sample_name")

ABI_metadata <- as.data.frame(ABI_metadata)
ABI_data <- t(as.matrix(ABI_data))
ABI_data <- as(ABI_data, "dgCMatrix")

print("Creating Seurat")
ABI_data <- CreateSeuratObject(ABI_data, meta.data = ABI_metadata)

print("Normalizing ABI")
ABI_data <- NormalizeData(ABI_data, normalization.method = "LogNormalize", scale.factor = 10000)

print("Finding variable genes")
ABI_data <- FindVariableFeatures(ABI_data, selection.method = "vst", nfeatures = 3000)
ABI_var_genes <- VariableFeatures(ABI_data)

print("Creating SCE")
ABI_data@meta.data["study_id"] <- rep("ABI", dim(ABI_data)[2])
ABI_data <- as.SingleCellExperiment(ABI_data)
ABI_data <- ABI_data[ABI_var_genes,]

print("ABI SCE created")

snRNA <- readRDS("/home/anjali5/projects/def-cnagy/1_enhanced_harmonized_object.Rds")
snRNA@meta.data["study_id"] <- rep("snRNA", dim(snRNA)[2])
snRNA
#snRNA@meta.data['Cluster_combined'] <- plyr::mapvalues(snRNA@meta.data$Cluster, from = levels(snRNA@meta.data$Cluster), to = c("InPV", "InPV", "ExN", "ExN_L56", "ExN_L56", "ExN_L56",  "ExN_L56", "ExN_L56",  "ExN", "ExN","ExN_L24","ExN_L46", "ExN","ExN_L23" , "ExN_L23", "ExN_L24", "ExN_L46",  "Ast1", "Ast2", "InN_ADARB2", "InN_ADARB2", "InN_VIP", "InN_VIP", "ExN_L56", "ExN_L56",  "InN_SST", "InN_SST", "InN_LAMP5", "Mix", "OPC3", "Oli3","Oli1", "Oli2", "OPC1", "OPC2", "ExN", "ExN_L35", "ExN", "Mix", "Endo", "Micro"))
snRNA <- subset(snRNA, subset=Broad!="Mix")
snRNA <- FindVariableFeatures(snRNA, selection.method = "vst", nfeatures = 3000)
snRNA_var_genes <- VariableFeatures(snRNA)

snRNA <- DietSeurat(snRNA, graphs = "pca")
snRNA <- as.SingleCellExperiment(snRNA)
snRNA <- snRNA[snRNA_var_genes,]

#snRNA
#print(snRNA@assays)

print("Training snRNA on Sub cluster")

pretrained_model_ABI = MetaNeighbor::trainModel(
 var_genes = ABI_var_genes,
 dat = ABI_data,
 study_id = ABI_data$study_id,
 cell_type = ABI_data$subclass_label
)


pretrained_model_snRNA = MetaNeighbor::trainModel(
 var_genes = snRNA_var_genes,
 dat = snRNA,
 study_id = snRNA$study_id,
 cell_type = snRNA$Broad
)

#save model to file

write.table(pretrained_model_snRNA,
        "/home/anjali5/scratch/MetaNeighbor/pretrained_snRNA_broad_promoter")
write.table(pretrained_model_ABI,
        "/home/anjali5/scratch/MetaNeighbor/pretrained_ABI_broad_promoter")

print("saved")
rm(ABI_data)
rm(snRNA)

print("snATC data")
ArchR@meta.data["study_id"] <- rep("snATAC", dim(ArchR)[2])

print("Convert to SCE")
ArchR <- as.SingleCellExperiment(ArchR,assay = "RNA")
#save.image("~/scratch/metaneighbor_tmp.rds")

plots_list <- list()

print("Running comparisons")
pdf(file = "/home/anjali5/scratch/MetaNeighbor/snRNA_ABI_heatmaps_broad_3000_promoter.pdf",
        onefile = TRUE)
for(model_name in c("pretrained_model_snRNA", "pretrained_model_ABI")) {
	print(model_name)
	#aurocs
	aurocs = MetaNeighborUS(
	trained_model = get(model_name), 
	dat = ArchR,
	study_id = ArchR$study_id,
	cell_type = colData(ArchR)["ClustersMapped"][[1]],
	fast_version = TRUE)

	plotHeatmapPretrained(aurocs, cex = 0.7)

	#best hits
	best_hits = MetaNeighborUS(
	trained_model = get(model_name),
	dat = ArchR,
	study_id = ArchR$study_id,
	cell_type = colData(ArchR)["ClustersMapped"][[1]],
		one_vs_best = TRUE,symmetric_output = FALSE,
	        fast_version = TRUE)
	
	plotHeatmapPretrained(best_hits, cex = 1.5)

	#find and save top hits
	top_hits = topHitsByStudy(aurocs,
                threshold = 0)

	write.table(top_hits, paste0("/home/anjali5/scratch/MetaNeighbor/",
	  model_name,"_top_hits_for_ArchR_broad.csv"))
}

dev.off()

print(summary(warnings()))
sessionInfo()

#save.image("~/scratch/metaneighbor.rds")
#print("saved")
