##Clustering QC

library(ArchR)
library(Signac)
library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

##Load ArchR Project
proj <- loadArchRProject(path = "/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded")

######IMPORTNAT###### #Order your project alphabetically for table####
newdata <- proj[order(proj$Subject),]


#PCA of cluster proportions for each subject coloured by Batch, Sex, Condition
plots_list <- list()
subject_cluster_proportions <- table(newdata$Subject, newdata$Clusters)/rowSums(table(newdata$Subject, newdata$Clusters))
data_for_pca <- data.frame(Subject = unique(newdata$Subject),
                           Batch = stringr:: str_split(unique(paste(newdata$Subject, newdata$batch)), pattern = " ", simplify = TRUE)[,2],
                           Sex =  stringr:: str_split(unique(paste(newdata$Subject, newdata$sex)), pattern = " ", simplify = TRUE)[,2],
                           Condition =  stringr:: str_split(unique(paste(newdata$Subject, newdata$condition)), pattern = " ", simplify = TRUE)[,2])
res_pca <- prcomp(subject_cluster_proportions, scale. = TRUE) 
colour_vars <- c("Batch", "Sex","Condition")
for(colour_var in colour_vars) {
  plots_list[[colour_var]] <- ggbiplot::ggbiplot(res_pca, groups = data_for_pca[,colour_var], ellipse = TRUE, var.axes = FALSE) + geom_text(label=data_for_pca[,"Subject"])
}
pdf(file = paste0("/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/Subject_cluster_proportions_PCA_WithoutBatchCorrection_.pdf"),
    onefile = TRUE, height = 10, width = 12)
plots_list
dev.off()

##Heatmap for variable proportions in each cluster
confounders <- c("Subject", "batch", "condition", "sex") 
for(var in confounders) {
  tmp <- getCellColData(newdata, select = var, drop = TRUE)
  frequencies <- table(tmp, newdata$Clusters)
  percent_pre <- colSums(frequencies>0)/dim(frequencies)[1]
  #if(var %in% c("nFrags", "TSSEnrichment")){
  #median_per_cluster <- aggregate(proj$var, list(proj$Clusters), FUN=median)
  pdf(file = paste0("/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/Heatmap_confounders_Retained_Clusters_", var, "_.pdf"),
      onefile = TRUE, height = 13, width = 13)
  pheatmap::pheatmap(frequencies[,colnames(frequencies)[order(as.numeric(colnames(frequencies)))]], scale = "row")	
}
dev.off()


#PCA of subject proportions for each cluster
colour_vars <- c("Subject","batch", "sex","condition")
plots_list <- list()
for(colour_var in colour_vars) {
  tmp <- getCellColData(newdata, select = colour_var, drop = TRUE)
  cluster_proportions <- table(newdata$Clusters, tmp)/rowSums(table(newdata$Clusters, tmp))
  res_pca <- prcomp(na.omit(cluster_proportions), scale. = TRUE)
  data_for_pca <- cbind(as.data.frame(res_pca$x), clusters = rownames(res_pca$x))
  plots_list[[colour_var]] <- ggplot(data = data_for_pca, aes(x=PC1, y = PC2, label = clusters)) + geom_text() 
}
pdf(file = paste0("/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/Cluster_subject_proportions_PCA_.pdf"),
    onefile = TRUE, height = 10, width = 12)
plots_list
dev.off()
#Remove Outlier cluster & Repeat PCA for remaining clusters


#Plot basic QC parameters per Sample & Cluster
plots_list_cluster <- list()
plots_list_sample <- list()
features_to_plot <- c("nFrags", "TSSEnrichment", "mito_percent", "BlacklistRatio", "DoubletEnrichment", "amulet")

for(feature in features_to_plot) {
  plots_list_cluster[[paste0("Clusters_", feature)]] <- plotGroups(newdata, groupBy = "Clusters", name = feature, log2Norm = TRUE, plotAs = "violin", maxCells = 5000)
  plots_list_sample[[paste0("Sample_", feature)]] <- plotGroups(newdata, groupBy = "Sample", name = feature, log2Norm = TRUE, plotAs = "violin", maxCells = 5000)
}
pdf(file = "/home/anjali5/projects/def-cnagy/anjali5/Finalized_outputs/Per_Sample_Cluster_nFrags_TSS_mito_blacklist_doublets_.pdf",
    onefile = TRUE, height = 13, width = 13)
plots_list_cluster
plots_list_sample
dev.off()

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

##Ratios of Case-Control & Male-Female by Clusters
print("ratios case/control")

frequencies <- data.frame(table(proj$condition, proj$Cluster))
Ratios <- frequencies %>% group_by(Var2) %>% mutate(ratio = Freq[1]/Freq[2])
unique(Ratios$Var2)
unique(Ratios$ratio)

print("ratios male/female")
frequencies <- data.frame(table(proj$sex, proj$Cluster))
Ratios <- frequencies %>% group_by(Var2) %>% mutate(ratio = Freq[2]/Freq[1])
unique(Ratios$ratio)

df <- data.frame(proj@cellColData)
df_male <- dplyr::filter(df, sex=="male")
df_female <- dplyr::filter(df, sex=="female")

print("ratios male case/control")
frequencies <- data.frame(table(df_male$condition, df_male$Cluster))
Ratios <- frequencies %>% group_by(Var2) %>% mutate(ratio = Freq[1]/Freq[2])
unique(Ratios$ratio)

print("ratios female case/control")
frequencies <- data.frame(table(df_female$condition, df_female$Cluster))
Ratios <- frequencies %>% group_by(Var2) %>% mutate(ratio = Freq[1]/Freq[2])
unique(Ratios$ratio)

print("Done..")

#####Correlate variables (Sample, batch, sex) with UMAP reduced dimesions###
df <- proj@cellColData
rd <- getReducedDims(proj) #iterative LSI (before batch correction)
rd2 <- getReducedDims(proj, reducedDims = "Harmony") #(after batch correction)
df <- df[rownames(df) %in% rownames(rd),]
dim(df)
for(var in c("Subject", "batch", "sex")){
  print(var)
  UMAP_resid1 <- lm(rd[,1] ~ get(var), data=df, na.action=na.omit)
  UMAP_resid2 <- lm(rd[,2] ~ get(var), data=df, na.action=na.omit)
  Harmony_resid1 <- lm(rd2[,1] ~ get(var), data=df, na.action=na.omit)
  Harmony_resid2 <- lm(rd2[,2] ~ get(var), data=df, na.action=na.omit)
  print(summary(UMAP_resid1)$r.squared)
  print(summary(UMAP_resid2)$r.squared)
  print(summary(Harmony_resid1)$r.squared)
  print(summary(Harmony_resid2)$r.squared)
}

print(summary(warnings()))
sessionInfo()






