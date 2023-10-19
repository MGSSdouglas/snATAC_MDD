
#MetaNeighbor

library(MetaNeighbor)
library(Seurat)
library(SingleCellExperiment)
library(ArchR)

library(tidyverse)
sessionInfo()
traceback()

print("Loading data")
ArchR <- readRDS("~/projects/def-gturecki/anjali5/ArchR_Subcluster_by_peaks_updated.rds")
ArchR$SubClusters <- plyr::mapvalues(ArchR$SubClusters, from = unique(ArchR$SubClusters), to = c("ExN12","ExN9", "Mix2","ExN2", "ExN6", "ExN10","ExN8", "ExN11", "ExN1","InPV", "Oli7", "Ast2", "Oli5","ExN4", "ExN5","ExN3", "InLAMP5","InSST", "OPC2", "InVIP", "Oli1", "Oli4", "End1", "Ast3", "OPC4", "End2", "Mic1", "Mic2", "Oli3", "Mix1", "Oli2", "Ast4", "ExN7", "OPC1", "OPC3", "InN3", "OLi6", "Ast1"))

#load Jakel snRNA-seq data
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118257

jakel <- readRDS("/home/anjali5/projects/def-gturecki/anjali5/jakel.rds")
jakel <- FindVariableFeatures(jakel, nfeatures = 3000)
jakel_var_genes <- VariableFeatures(jakel)

print("Creating SCE")
jakel@meta.data["study_id"] <- rep("J", dim(jakel)[2])

jakel <- DietSeurat(jakel, graphs = "pca")
jakel <- as.SingleCellExperiment(jakel)
jakel <- jakel[jakel_var_genes,]

print("Jakel SCE created")

print("Training")

pretrained_model_jakel = MetaNeighbor::trainModel(
 var_genes = jakel_var_genes,
 dat = jakel,
 study_id = jakel$study_id,
 cell_type = jakel$celltype
)

#save model to file

write.table(pretrained_model_jakel,
        "/home/anjali5/scratch/MetaNeighbor/pretrained_jakel_celltype_3000")

print("snATC data")
ArchR@meta.data["study_id"] <- rep("snATAC", dim(ArchR)[2])


print("Convert to SCE")
ArchR <- as.SingleCellExperiment(ArchR,assay = "RNA")

plots_list <- list()

print("Running comparisons")
pdf(file = "/home/anjali5/scratch/MetaNeighbor/Jakel_snRNA_heatmaps_SubClusters_Promoter_3000.pdf",
        onefile = TRUE)
for(model_name in c("pretrained_model_jakel")) {
	print(model_name)
	#aurocs
	aurocs = MetaNeighborUS(
	trained_model = get(model_name), 
	dat = ArchR,
	study_id = ArchR$study_id,
	cell_type = colData(ArchR)["SubClusters"][[1]],
	fast_version = TRUE)

	plotHeatmapPretrained(aurocs, cex = 0.7)

	#best hits
	best_hits = MetaNeighborUS(
	trained_model = get(model_name),
	dat = ArchR,
	study_id = ArchR$study_id,
	cell_type = colData(ArchR)["SubClusters"][[1]],
		one_vs_best = TRUE,symmetric_output = FALSE,
	        fast_version = TRUE)
	
	plotHeatmapPretrained(best_hits, cex = 0.7)

	#find and save top hits
	top_hits = topHitsByStudy(aurocs,
                threshold = 0)

	write.table(top_hits, paste0("/home/anjali5/scratch/MetaNeighbor/",
	  model_name,"_top_hits_for_ArchR_SubClusters_promoter_3000.csv"))
}

dev.off()

print(summary(warnings()))
sessionInfo()


