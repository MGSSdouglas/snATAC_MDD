library(GenomicRanges)

#synapse: https://www.synapse.org/#!Synapse:syn12245061/tables/
#histone modification data - read and reduce

filenames <- list.files("~/Downloads/Neuronal_H3K4me3/", pattern="*.bed", full.names=TRUE)
ldf <- lapply(filenames, function(x) read.table(x,colClasses = c(rep("character", 3), rep("NULL", 7)),header=FALSE,col.names=c("chr","start","end",rep("na",7))))
res <- lapply(ldf, function(x) GenomicRanges::makeGRangesFromDataFrame(x))
combined.peaks.neuron.h3k4me3 <- GenomicRanges::reduce(x = c(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],res[[6]],res[[7]],res[[8]],res[[9]],res[[10]],res[[11]],res[[12]],res[[13]],res[[14]],res[[15]],res[[16]],res[[17]]))
write.table(combined.peaks.neuron.h3k4me3, "~/Desktop/NeuN+_H3K4me3.bed", col.names = F, row.names = F, quote=F)

filenames <- list.files("~/Downloads/Neuronal_H3K27ac/", pattern="*.bed", full.names=TRUE)
ldf <- lapply(filenames, function(x) read.table(x,colClasses = c(rep("character", 3), rep("NULL", 7)),header=FALSE,col.names=c("chr","start","end",rep("na",7))))
res <- lapply(ldf, function(x) GenomicRanges::makeGRangesFromDataFrame(x))
combined.peaks.neuron.h3k4me3 <- GenomicRanges::reduce(x = c(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],res[[6]],res[[7]],res[[8]],res[[9]],res[[10]],res[[11]],res[[12]],res[[13]],res[[14]],res[[15]],res[[16]],res[[17]]))
write.table(combined.peaks.neuron.h3k4me3, "~/Desktop/NeuN+_H3K27ac.bed", col.names = F, row.names = F, quote=F)


filenames <- list.files("~/Downloads/Glia_H3K27ac/", pattern="*.bed", full.names=TRUE)
ldf <- lapply(filenames, function(x) read.table(x,colClasses = c(rep("character", 3), rep("NULL", 7)),header=FALSE,col.names=c("chr","start","end",rep("na",7))))
res <- lapply(ldf, function(x) GenomicRanges::makeGRangesFromDataFrame(x))
combined.peaks.glia.h3k27ac <- GenomicRanges::reduce(x = c(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],res[[6]],res[[7]],res[[8]],res[[9]],res[[10]],res[[11]],res[[12]],res[[13]],res[[14]],res[[15]],res[[16]],res[[17]]))
write.table(combined.peaks.glia.h3k27ac, "~/Desktop/NeuN-_H3K27ac.bed", col.names = F, row.names = F, quote=F)

filenames <- list.files("~/Downloads/Glia_H3K4me3/", pattern="*.bed", full.names=TRUE)
ldf <- lapply(filenames, function(x) read.table(x,colClasses = c(rep("character", 3), rep("NULL", 7)),header=FALSE,col.names=c("chr","start","end",rep("na",7))))
res <- lapply(ldf, function(x) GenomicRanges::makeGRangesFromDataFrame(x))
combined.peaks.glia.h3k4me3 <- GenomicRanges::reduce(x = c(res[[1]],res[[2]],res[[3]],res[[4]],res[[5]],res[[6]],res[[7]],res[[8]],res[[9]],res[[10]],res[[11]],res[[12]],res[[13]],res[[14]],res[[15]],res[[16]],res[[17]]))
write.table(combined.peaks.glia.h3k4me3, "~/Desktop/NeuN-_H3K4me3.bed", col.names = F, row.names = F, quote=F)

