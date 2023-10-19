# source("3__ensid_to_gene_symbol.R")

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
processed_data_path <- paste(GKMSVM_PREPARED_DATA_DIR, data_id, sep = "/")
dir.create(processed_data_path, showWarnings = FALSE, recursive = TRUE)

# ensembl details 
ensembl <- useEnsembl(biomart = "genes", host = "https://useast.ensembl.org")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

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
