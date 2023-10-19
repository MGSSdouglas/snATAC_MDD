# source("2__calculate_FDR_and_subset.R")

library(data.table)
library(biomaRt)

###### configuration

# get environment variable values 
GKMSVM_WORKSPACE_DIR <- Sys.getenv("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR <- Sys.getenv("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR <- Sys.getenv("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR <- Sys.getenv("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR <- Sys.getenv("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR <- Sys.getenv("GKMSVM_BIN_DIR")

# configure paths
data_id <- "eqtl_catalogue"
raw_data_dir <- paste(GKMSVM_RAW_DATA_DIR, data_id, sep = "/")
sumstats_dir <- paste(raw_data_dir, "sumstats", sep = "/")
processed_data_path <- paste(GKMSVM_PREPARED_DATA_DIR, data_id, sep = "/")
dir.create(processed_data_path, showWarnings = FALSE, recursive = TRUE)

###### load sumstats, calculate FDR, subset to ssnp and cdsnps and save

### BrainSeq - gene expression

# read sumstats
filepath <- "QTD000051.all.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "BrainSeq.ge.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.BrainSeq.ge.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.BrainSeq.ge.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.BrainSeq.ge.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished BrainSeq - gene expression")


### BrainSeq - splicing

# read sumstats
filepath <- "QTD000055.cc.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "BrainSeq.sp.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")


print("Finished BrainSeq - splicing")


### CommonMind - gene expression

# read sumstats
filepath <- "QTD000075.all.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "CommonMind.ge.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished CommonMind - gene expression")


### CommonMind - splicing

# read sumstats
filepath <- "QTD000079.cc.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "CommonMind.sp.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished CommonMind - splicing")


### GTEx - gene expression

# read sumstats
filepath <- "QTD000176.all.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "GTEx.ge.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished GTEx - gene expression")


### GTEx - splicing

# read sumstats
filepath <- "QTD000180.cc.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "GTEx.sp.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")


print("Finished GTEx - splicing")


### ROSMAP - gene expression

# read sumstats
filepath <- "QTD000434.all.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "ROSMAP.ge.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished ROSMAP - gene expression")


### ROSMAP - splicing

# read sumstats
filepath <- "QTD000438.cc.tsv.gz"
filepath <- paste(sumstats_dir, filepath, sep = "/")
sumstats <- fread(filepath)

# calculate FDR and retain FDR < 0.05
sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
sumstats <- subset(sumstats, fdr < 0.05)

# save FDR thresholded sumstats
filepath <- paste(processed_data_path, "ROSMAP.sp.tsv", sep = "/")
fwrite(sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished ROSMAP - splicing")


print("Script finished")
