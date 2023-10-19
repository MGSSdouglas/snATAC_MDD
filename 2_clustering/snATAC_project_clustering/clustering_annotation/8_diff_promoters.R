library(Seurat)
library(Signac)

brain <- readRDS("~/projects/def-gturecki/anjali5/ArchR_cluster_by_peaks_updated.rds")
#atac
#Idents(atac) <- "ClustersMapped"
#brain <- subset(x=atac, idents=c("Mic1", "Mic2", "End1"))
#brain <- subset(x=atac, idents=c("Ast1","Ast2","Ast3","Ast4"))
#brain <- subset(x=atac, idents="ExN")
brain
#print("chromVar saved")
DefaultAssay(brain) <- 'RNA'

diff <- data.frame()
print("running diff..")
Idents(brain) <- "Clusters"
#loop
differential.activity <- FindAllMarkers(
  object = brain,assay="RNA",test.use = 'LR',
  min.pct = 0.1,only.pos = TRUE,
  latent.vars = "nCount_RNA"
)
write.csv(differential.activity,"~/scratch/All_clusters_diffmarkers.csv")
