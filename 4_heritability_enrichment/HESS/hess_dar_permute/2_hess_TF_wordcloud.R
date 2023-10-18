library(wordcloud2)
library(wordcloud)

#TF motifs enriched in marker peaks in significant LD blocks
marker <- readxl::read_excel("~/Ex_marker_hess.xlsx")

#TF motifs enriched in DARs in significant LD blocks
dar <- read.csv("~/ExN1_DAR_knownResults.txt",sep="\t")
dar <- dar[dar$q.value..Benjamini. < 0.05,]
marker <- marker[marker$`q-value (Benjamini)` < 0.05,]

#get TF families common to both
common <- intersect(marker$`Motif Name`, dar$Motif.Name)
tf <- do.call(rbind,stringr::str_split(common, "/"))[,1]
families <- gsub("[\\(\\)]", "", regmatches(tf, gregexpr("\\(.*?\\)", tf)))
families[41] <- "bHLH"
df <- data.frame(table(families))

#plot TF families based on frequency of individual TF members enriched
wordcloud(words = df$families, freq = df$Freq, min.freq = 1,max.words=200, random.order=FALSE, rot.per=0, colors=brewer.pal(8, "Dark2"))
