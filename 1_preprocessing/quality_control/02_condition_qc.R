library(ArchR)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggplot2)

proj <- loadArchRProject("~/projects/def-gturecki/ArchR_output/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered")

proj@cellColData %>% as.data.frame() %>% group_by(Subject) %>% summarise(nFrags=mean(nFrags), TSSEnrichment=mean(TSSEnrichment),
                                                     pH=mean(as.numeric(pHValue)), PMI=mean(as.numeric(PMI)), Age=mean(as.numeric(Age)), FRIP=mean(as.numeric(FRIP)), condition=unique(condition)) -> df

features_to_plot <- c("nFrags", "TSSEnrichment", "FRIP", "Age", "PMI", "pH")
plot_list <- list()
for(feature in features_to_plot){
  print(feature)
  #values=df[[feature]]
  p <- ggplot(df, aes(x=condition,y=.data[[feature]],fill=condition)) + geom_violin() + stat_compare_means(method="wilcox.test", label.x.npc = "center") +theme_classic() + scale_fill_manual(values=c("salmon4","tan")) + theme(text = element_text(size = 20),axis.text = element_text(size = 20))+ ylab(feature) + geom_jitter(position = position_jitter(seed = 1)) 
  plot_list[[feature]] <- p
  #plots_list[[feature]] <- ggplot(df, aes(x=condition,y=df[[feature]],fill=condition, label=Subject)) + geom_violin() + stat_compare_means(method="wilcox.test", label.x.npc = "center") +theme_classic() + scale_fill_manual(values=c("salmon4","tan")) + theme(text = element_text(size = 20),axis.text = element_text(size = 20)) + geom_text(check_overlap = TRUE, position=position_jitter(width=0.09),label.size = 0.1) + ylab(feature)
}

pdf(file = "/home/anjali5/projects/def-gturecki/anjali5/Per_Condition_QC_stats.pdf",
    onefile = TRUE)
plot_list
dev.off()

#get mean and SE
result <- df %>%
  group_by(condition) %>%
  summarize(
    Mean_age = round(mean(Age), 2),
    SE_age = sd(Age) / sqrt(n()),  
    Mean_pmi = round(mean(PMI),2),
    SE_pmi = sd(PMI) / sqrt(n()),
    Mean_ph = round(mean(pH),2),
    SE_ph = sd(pH) / sqrt(n())
  )
print(result)
# A tibble: 2 x 7
#condition Mean_age SE_age Mean_pmi SE_pmi Mean_ph  SE_ph
#<chr>        <dbl>  <dbl>    <dbl>  <dbl>   <dbl>  <dbl>
#  1 Case          44.6   2.51     44.2   2.70    6.54 0.0471
#2 Control       44.2   2.75     33.4   3.57    6.43 0.0408


