library(ArchR)

proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered")
proj

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Subjects", 
  geneSymbol = "XIST", 
  upstream = 1000,
  downstream = 1000
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-per-subject-XIST.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 15, height = 30)

#Plot sex-specific genes per subject
plots_list <- list()
features_to_plot <- c("XIST", "NLGN4Y", "JPX", "UTY")
for(feature in features_to_plot) {
  plots_list[[feature]] <- plotGroups(newdata, groupBy = "Subject",colorBy="GeneScoreMatrix", name = feature, log2Norm = TRUE, plotAs = "violin", maxCells = 5000)
}
pdf(file = "/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/Per_Subject_Sex_Validation_genes_.pdf",
    onefile = TRUE, height = 13, width = 13)
plots_list
dev.off()
