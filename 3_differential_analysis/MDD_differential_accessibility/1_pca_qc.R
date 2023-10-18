library(PCAtools)
library(DESeq2)
#x1 <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
x1 <- paste0("PC",1:84)
plots_list <- list()

for(r in 2:length(x1)-1){
  print(r)
  plots_list[[paste0("pca",r)]] <- PCAtools::biplot(p,x=x1[r], y=x1[r+1],colby="condition",drawConnectors = FALSE,legendPosition = "bottom",lab = NULL)
}

pdf(file = "~/plot_PC_condition_84PC.pdf",onefile = TRUE)
plots_list
dev.off()

pb <- readRDS("~/pb_sample_aggr.rds")
obj <- readRDS("~/seuratObj_PeakMatrix_ArchR.rds") #, ReadsInPeaks = sum(as.numeric(ReadsInPeaks)
tibble(obj@meta.data) %>% group_by(Subject) %>% summarise(
  age = unique(as.numeric(Age)),condition=unique(condition),
  pH = unique(as.numeric(pHValue)),nFrags = mean(as.numeric(nFrags)),
  pmi = unique(as.numeric(PMI)), 
  batch = unique(batch), sex=unique(sex)) %>% ungroup() -> meta
obj@meta.data %>% group_by(Subject) %>% summarise(count=n()) -> counts_df
meta$nCell <- counts_df$count
meta <- data.frame(meta)
rownames(meta) <- meta$Subject
ncol(pb)==nrow(meta)
dds <- DESeqDataSetFromMatrix(countData=assay(pb),colData=meta, tidy = FALSE, design = ~batch+age+sex+pmi+nCell+condition+nFrags)
#Normalise the data and transform the normalised counts to variance-stabilised expression levels:
dds <- DESeq(dds)
vst <- assay(vst(dds))
p <- pca(vst, metadata = colData(dds), removeVar = 0)
screeplot(p,drawCumulativeSumLine = FALSE,components = 10)
eigencorplot(p,components=getComponents(p, seq_len(20)),metavars = c("condition", "age", "sex","batch","nCell","pmi","pH"))


my_list=list()
for(i in 1:42){
  d[[i]]@path -> f
  my_list[[length(my_list) + 1]] <- f
}


