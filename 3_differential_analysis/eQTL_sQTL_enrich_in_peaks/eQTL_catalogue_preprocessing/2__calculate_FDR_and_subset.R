library(data.table)
library(biomaRt)

### This script loads the downloaded summary statistics, calculates FDR and retains (SNP, gene) pairs having FDR < 0.05.

###### configuration
raw_data_dir <- "./data/raw"
sumstats_dir <- "./data/raw/sumstats"

processed_data_path <- "./data/processed"
dir.create(processed_data_path, showWarnings = FALSE, recursive = TRUE)

###### load MDD relevant, cdsnps and ssnps 

# get unique MDD relevant SNPs
filepath <- "MDD_relevant_snps.tsv"
filepath <- paste(raw_data_dir, filepath, sep = "/")
mdd_snps <- fread(filepath, header = T, sep = "\t")
mdd_snps <- unique(mdd_snps$Name)

# get unique cdSNPs
filepath <- "cdSNP_table.no_interaction.final.51bp_null.with_normalization.v13.tsv"
filepath <- paste(raw_data_dir, filepath, sep = "/")
cdsnps <- fread(filepath, header = T, sep = "\t")
cdsnp_rsids <- unique(cdsnps$snp__Name)

# get unique sSNPs
filepath <- "sSNP_table.no_interaction.final.51bp_null.with_normalization.v13.tsv"
filepath <- paste(raw_data_dir, filepath, sep = "/")
ssnp <- fread(filepath, header = T, sep = "\t")
ssnp_rsids <- unique(ssnp$snp__Name)

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.BrainSeq.sp.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.BrainSeq.sp.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.BrainSeq.sp.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.CommonMind.ge.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.CommonMind.ge.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.CommonMind.ge.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.CommonMind.sp.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.CommonMind.sp.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.CommonMind.sp.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.GTEx.ge.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.GTEx.ge.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.GTEx.ge.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.GTEx.sp.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.GTEx.sp.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.GTEx.sp.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.ROSMAP.ge.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.ROSMAP.ge.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.ROSMAP.ge.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

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

# retain ssnps in sumstats and save
ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
filepath <- paste(processed_data_path, "sSNP.ROSMAP.sp.tsv", sep = "/")
fwrite(ssnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain cdsnps in sumstats and save
cdsnp_sumstats <- subset(sumstats, rsid %in% cdsnp_rsids)
filepath <- paste(processed_data_path, "cdSNP.ROSMAP.sp.tsv", sep = "/")
fwrite(cdsnp_sumstats, filepath, sep = "\t", quote = F, na = "NA")

# retain MDD relevant SNPs in sumstats and save
mdd_sumstats <- subset(sumstats, rsid %in% mdd_snps)
filepath <- paste(processed_data_path, "MDD.ROSMAP.sp.tsv", sep = "/")
fwrite(mdd_sumstats, filepath, sep = "\t", quote = F, na = "NA")

print("Finished ROSMAP - splicing")


# ### GTEx v8 imported data - gene expression

# # read sumstats
# filepath <- "Brain_Frontal_Cortex_BA9.tsv.gz"
# filepath <- paste(sumstats_dir, filepath, sep = "/")
# sumstats <- fread(filepath)

# # calculate FDR and retain FDR < 0.05
# sumstats$fdr <- p.adjust(sumstats$pvalue, 'BH')
# sumstats <- subset(sumstats, fdr < 0.05)

# # retain ssnps in sumstats and save
# ssnp_sumstats <- subset(sumstats, rsid %in% ssnp_rsids)
# filepath <- paste(processed_data_path, "sSNP.GTExv8.ge.tsv", sep = "/")

print("Script finished")
