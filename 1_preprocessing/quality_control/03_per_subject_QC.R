library(ArchR)

proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")
proj

###compute per subject qc

p1 <- plotFragmentSizes(ArchRProj = proj,groupBy = "Subject",threads =32)

saveRDS(p1, "~/scratch/FigfragmetSizeDistribution.rds")

print("done")

p2 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Subject",
    colorBy = "cellColData",
    name = "nFrags",
    plotAs = "ridges"
   )

saveRDS(p2, "~/scratch/FignFrags.rds")

p3 <- plotTSSEnrichment(ArchRProj = proj,groupBy = "Subject",threads =32)

saveRDS(p3, "~/scratch/FigTSSEnrichemnt.rds")


df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
saveRDS(p, "~/scratch/AllQC.rds")

p <- readRDS("~/scratch/FigfragmetSizeDistribution.rds")
pdf("~/projects/def-gturecki/anjali5/SuppFigData/FigfragmetSizeDistribution.pdf")
p
dev.off()

p <- readRDS("~/scratch/FigTSSEnrichemnt.rds")
pdf("~/projects/def-gturecki/anjali5/SuppFigData/FigTSSEnrichemnt.pdf")
p
dev.off()

p <- readRDS("~/scratch/FignFrags.rds")
pdf("~/projects/def-gturecki/anjali5/SuppFigData/FignFrags.pdf")
p
dev.off()

p <- readRDS("~/scratch/AllQC.rds")
pdf("~/projects/def-gturecki/anjali5/SuppFigData/AllQC.pdf")
p
dev.off()
#plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Subject",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )


saveRDS(p2, "~/scratch/FigTSS.rds")
p <- readRDS("~/scratch/FigTSS.rds")
pdf("~/projects/def-gturecki/anjali5/SuppFigData/FigTSS.pdf")
p
dev.off()
