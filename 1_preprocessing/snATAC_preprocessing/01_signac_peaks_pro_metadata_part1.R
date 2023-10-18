library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(GenomeInfoDb)
library(data.table)

#####Call bash args########
args <- commandArgs(TRUE)

subject <- noquote(args[1])

print(args)
print(subject)

Dir="/home/anjali5/scratch/graham_atac_analysis/"

#####fragmets#####
frags <- CreateFragmentObject(
path = paste0(Dir,subject,"/summary/fragments.tsv.gz"),
cells = NULL
)
#####Call_Peaks#####
peaks <- CallPeaks(
  object = frags,
  macs2.path = "/home/anjali5/.local/bin/macs2",cleanup = TRUE
)
print(head(peaks))
print(length(peaks))
#####Sanitize peaks####
chr_sizes_file="/home/anjali5/scATAC-pro/annotation/chrom_hg38.sizes"
chrs = fread(chr_sizes_file, header = F)$V1
peaks <- peaks[peaks@seqnames %in% chrs]
peaks <- subsetByOverlaps(peaks,blacklist_hg38_unified,invert = TRUE)
print(head(peaks))
print(length(peaks))
##SAVE############ 
peak_file=paste0(Dir,subject,"/peaks/MACS2/peaks.bed")
write.table(peaks, file = peak_file, sep = '\t', row.names = F, col.names = F, quote = F)
print("peaks done")
