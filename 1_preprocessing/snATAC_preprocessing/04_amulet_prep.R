library(dplyr)
library(Matrix)
library(data.table)

####pass argumets to R##
args <- commandArgs(TRUE)
subject <- noquote(args[1])
print(args)

Dir="/home/anjali5/scratch/graham_atac_analysis/"

###load and add to singlecell.csv
md10x <- read.csv(
file = paste0(Dir,"atac_",subject,"/singlecell.csv"),
sep = ",",
header = TRUE
)

###add cellranger_cell for 10x-identified cells and is__cell_barcode=1 for amulet
barcodes <- read.csv(paste0(Dir,"atac_",subject,"/archR_barcodes"),header=TRUE,sep="")
md10x$cellranger_cell <- md10x$is__cell_barcode
print(paste0("saved is__cell_barcode as cellranger_cell for:", subject))
md10x$is__cell_barcode <- ifelse(md10x$barcode %in% barcodes$barcode, "1", "0")
#md10x$cellranger_cell <- md10x$is__cell_barcode
#md10x$is__cell_barcode <- ifelse(md10x$is__cell_barcode=="0", "1", "1")
write.csv(md10x,paste0(Dir,"atac_",subject,"/singlecell_archR.csv"), quote=FALSE, sep=",",row.names=FALSE)
print("done")


