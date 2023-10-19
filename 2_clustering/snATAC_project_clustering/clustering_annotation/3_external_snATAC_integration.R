library(Matrix)
library(Signac)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(ArchR)
set.seed(1234)

#load Morabito et.al (reference) and snATAC data (query)
#https://www.synapse.org/#!Synapse:syn22079621/.
#ref
mbATAC <- readRDS("~/projects/def-gturecki/anjali5/Morabito_Validation/morabito_seuratObj.rds")
#query
obj <- readRDS("~/projects/def-gturecki/anjali5/ArchR_Subcluster_by_peaks_updated.rds")

#compute LSI on ref 
DefaultAssay(mbATAC) <- "peaks"
mbATAC <- FindTopFeatures(mbATAC, min.cutoff = 10)
mbATAC <- RunTFIDF(mbATAC)
mbATAC <- RunSVD(mbATAC)
mbATAC

#Prepare snATAC data
DefaultAssay(obj) <- "peaks"
frags <- Fragments(obj)
frags
cells <- colnames(obj)
meta  <- obj@meta.data

#Quantify ref snATAC peaks in query snATAC-seq dataset
counts <- FeatureMatrix(
  fragments = frags,
  features = mbATAC@assays$peaks@ranges,
  cells = cells
)
head(counts)

#print("counted")
#create Seurat object 
atac.assay <- CreateChromatinAssay(
  counts = counts,
  fragments = frags,genome=NULL
)

atac <- CreateSeuratObject(counts = atac.assay, assay = "peaks",meta.data = meta)
print("archr obj created..")

# compute LSI on snATAC-seq
atac <- FindTopFeatures(atac, min.cutoff = 10)
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30)

# first add dataset-identifying features
atac$dataset <- "ATAC"
mbATAC$dataset <- "Morabito"

#save.image("~/scratch/graham_merge_analysis/merged_Rdata/morabito_ArchR.Rdata")
# compute UMAP and store the UMAP model
mbATAC <- RunUMAP(mbATAC, reduction = "lsi", dims = 2:30, return.model = TRUE)

#find transfer anchors
print("find transfer anchors..")

transfer.anchors <- FindTransferAnchors(
  reference = mbATAC,
  query = atac,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)
print("FindTransferAnchors done..")
#print(head(transfer.anchors))

#map query onto the reference dataset
print("map query onto the reference dataset")

atac <- MapQuery(
  anchorset = transfer.anchors,
  reference = mbATAC,
  query = atac,
  refdata = mbATAC$celltype,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

print("MapQuery done..")
saveRDS(atac,"~/scratch/graham_merge_analysis/merged_Rdata/ArchR_predicted_by_morabito.rds")
print(median(atac$predicted.id.score))
print(mean(atac$predicted.id.score))

##plot..
pdf("~/projects/def-gturecki/anjali5/Finalized_outputs/morabito_cluster_onto_ArchR.pdf")
p1 <- DimPlot(mbATAC, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p2 <- DimPlot(atac, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query")

p1 | p2

dev.off()
print("Plotted..")

