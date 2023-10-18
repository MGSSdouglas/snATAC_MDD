library(cowplot)
#scree plot
rd <- getReducedDims(proj, reducedDims = "IterativeLSI",dimsToUse = 1:50)
sd <- colSds(rd) #colVars var <- colVars(as.matrix(rd))
ggplot(data.frame(dims = 1:50, stdev = sd[1:50]))+geom_point(mapping = aes_string(x = 'dims', y = 'stdev'))+labs(x = "LSI",y="Standard Deviation") + theme_cowplot()
sdr <- data.frame(rep("LSI",70),var)
colnames(sdr) <- c("lsi","sd")
sdr$lsi <- paste0(sdr$lsi,rownames(sdr))
lowest_sd <- sdr[order(sd),]  ##get LSI corresponding to least stdev
##since LSI_51 explains least variance selected 1:50 components for clustering##
#cluster by sample
df <- proj@cellColData
df %>% group_by(Clusters,Sample) %>% count() %>%  group_by(Clusters) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>%   ggplot(aes(x=Clusters,y=percent, fill=Sample)) + geom_col() +   ggtitle("Percentage of subjects  per cluster") + theme_cowplot()
#####
add_sub_clusters <- function(input, dimsToUse, resolution, name){
  proj <- addClusters(input = proj,reducedDims = "Harmony",name="Clusters",method = "Seurat",resolution = 0.2, k.param=20, algorithm=3, dimsToUse = 1:20,sampleCells=20000,force=TRUE)
  if(proj$Clusters)
    
    for(col in colnames(rd)){
      LSI <- multinom(rd[[col]]~proj$Clusters_2)
      print(paste(col,summary(LSI)$r.squared))
    }
}
