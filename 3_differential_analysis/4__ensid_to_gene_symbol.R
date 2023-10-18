library(data.table)
library(biomaRt)

### This script converts unique gene ensids to gene symbols for each study and modality type

###### configuration

processed_data_path <- "./data/processed"

# ensembl details 
ensembl <- useEnsembl(biomart = "genes", host = "https://useast.ensembl.org")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)



# ###### cdSNP

# snp_type <- "cdSNP"


# ### BrainSeq - ge
# study_name <- "BrainSeq"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



# ### BrainSeq - sp
# study_name <- "BrainSeq"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### CommonMind - ge

# study_name <- "CommonMind"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### CommonMind - sp
# study_name <- "CommonMind"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### GTEx - ge
# study_name <- "GTEx"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



# ### GTEx - sp
# study_name <- "GTEx"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### ROSMAP - ge
# study_name <- "ROSMAP"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### ROSMAP - sp
# study_name <- "ROSMAP"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")




# ###### sSNP

# snp_type <- "sSNP"


# ### BrainSeq - ge
# study_name <- "BrainSeq"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



# ### BrainSeq - sp
# study_name <- "BrainSeq"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### CommonMind - ge

# study_name <- "CommonMind"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### CommonMind - sp
# study_name <- "CommonMind"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### GTEx - ge
# study_name <- "GTEx"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



# ### GTEx - sp
# study_name <- "GTEx"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### ROSMAP - ge
# study_name <- "ROSMAP"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### ROSMAP - sp
# study_name <- "ROSMAP"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")







# ###### MDD relevant SNPs

# snp_type <- "MDD"


# ### BrainSeq - ge
# study_name <- "BrainSeq"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



# ### BrainSeq - sp
# study_name <- "BrainSeq"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### CommonMind - ge

# study_name <- "CommonMind"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### CommonMind - sp
# study_name <- "CommonMind"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### GTEx - ge
# study_name <- "GTEx"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



# ### GTEx - sp
# study_name <- "GTEx"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### ROSMAP - ge
# study_name <- "ROSMAP"
# modality_type <- "ge"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


# ### ROSMAP - sp
# study_name <- "ROSMAP"
# modality_type <- "sp"

# # read data 
# filepath <- paste(snp_type, study_name, modality_type, "tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# eqtls <- fread(filepath, header = T, sep = "\t")

# # get eqtl gene ensids
# eqtl_gene_ensids = unique(eqtls$gene_id)
# eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# # query gene symbols for eqtl gene ensids
# eqtl_gene_mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
#     filters = 'ensembl_gene_id',
#     values = unlist(eqtl_gene_ensids),
#     mart = ensembl
# )

# # save gene symbols for eqtl gene ensids
# filepath <- paste(snp_type, study_name, modality_type, "gene_symbol_mapping.tsv", sep = ".")
# filepath <- paste(processed_data_path, filepath, sep = "/")
# fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")






###### FDR < 0.05


### BrainSeq - ge
study_name <- "BrainSeq"
modality_type <- "ge"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



### BrainSeq - sp
study_name <- "BrainSeq"
modality_type <- "sp"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")

### CommonMind - ge

study_name <- "CommonMind"
modality_type <- "ge"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


### CommonMind - sp
study_name <- "CommonMind"
modality_type <- "sp"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


### GTEx - ge
study_name <- "GTEx"
modality_type <- "ge"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")



### GTEx - sp
study_name <- "GTEx"
modality_type <- "sp"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


### ROSMAP - ge
study_name <- "ROSMAP"
modality_type <- "ge"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")


### ROSMAP - sp
study_name <- "ROSMAP"
modality_type <- "sp"

# read data 
filepath <- paste(study_name, modality_type, "gene_ids.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
eqtls <- fread(filepath, header = T, sep = "\t")

# get eqtl gene ensids
eqtl_gene_ensids = unique(eqtls$gene_id)
eqtl_gene_ensids <- as.list(eqtl_gene_ensids)

# query gene symbols for eqtl gene ensids
eqtl_gene_mapping <- getBM(
    attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_accession', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = unlist(eqtl_gene_ensids),
    mart = ensembl
)

# save gene symbols for eqtl gene ensids
filepath <- paste(study_name, modality_type, "gene_id_mapping.tsv", sep = ".")
filepath <- paste(processed_data_path, filepath, sep = "/")
fwrite(eqtl_gene_mapping, file=filepath, sep = "\t")

print("Script finished")
