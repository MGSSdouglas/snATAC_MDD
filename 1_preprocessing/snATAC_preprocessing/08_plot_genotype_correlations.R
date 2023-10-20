library(dplyr)
library(Matrix)
library(pheatmap)

#load genotype correlations
genotype_match <- read.table("~/Downloads/cor_results_v2.csv",sep=",",header = TRUE, row.names = "MGSS_id")
length(colnames(genotype_match))
length(rownames(genotype_match))

#get matrix
mat <- genotype_match[,3:136]
mat <- mat[,colnames(mat) %in% rownames(mat)]
mat <- mat[,order(colnames(mat))]
mat <- mat[order(rownames(mat)),]
mat[,!colnames(mat) %in% c("S142", "S104")] -> mat
mat[!rownames(mat) %in% c("S142", "S104"),] -> mat
ifelse(!rownames(mat) %in% c("atac_A2B2female","atac_A4B4female","atac_A8B8female","atac_A14B14female","atac_A17B17female","atac_A18B18female", "atac_A19B19female", "atac_A23B23female", "atac_A28B28female","atac_A29B29female", "atac_A30B30female", "atac_A32B32female", "atac_A39B39female")
      , sub('.', '', rownames(mat)), rownames(mat)) -> replace_vector
rownames(mat) <- replace_vector
colnames(mat) <- rownames(mat)
rownames(mat)==colnames(mat) #sanity check

#load metadata
meta <- readxl::read_excel("/Users/anjalichawla/Desktop/psychencode_data/metadata_subjects_updated.xlsx")
#meta <- meta[meta$Subject!="atac_A32B32female",]

#map brainID to subject
for (i in 1:length(rownames(mat))) {
  currentsubject = rownames(mat)[i]
  print(currentsubject)
  if(currentsubject %in% meta$BrainID){
    rownames(mat)[i] <- meta$Subject[meta$BrainID==currentsubject]
  }
  else{
    print("Not available?")
  }
}
colnames(mat) <- rownames(mat)
pdf("~/Desktop/snatac_manuscript/genotype_match.pdf", height=15, width=15)
pheatmap::pheatmap(mat,cluster_rows = F,cluster_cols = F)
dev.off()
