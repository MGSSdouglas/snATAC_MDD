library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(tidyr)

#set workdir
setwd("/home/frosig/scratch/WASP_allelic_imbalance/create_sh_snRNA_snATAC_CHT")

# Load your SNP data
snp_data <- read.table("peaks84_dup_row.txt") #load the fixed ones (having each star-end into a single row)

#snp as list
snp_data <- read.table("peaks84.txt") #load the fixed ones (having each star-end into a single row)

# Load the GTF file to extract gene annotations
gtf_file <- "Homo_sapiens.GRCh38.109.gtf"
gtf_data <- rtracklayer::import(gtf_file)

# Filter for exons only
exon_data <- subset(gtf_data, type == "exon")
seqlevels(exon_data) <- paste0("chr",seqlevels(exon_data))


# Split exons by gene_id to retain gene information
exon_gr <- GRanges(seqnames = seqnames(exon_data), ranges = IRanges(start = start(exon_data), end = end(exon_data)), gene_id = exon_data$gene_id)
exons_by_gene <- split(exon_gr, exon_gr$gene_id)


# Collapse overlapping or duplicate exons within each gene
collapsed_exons_list <- lapply(exons_by_gene, function(x) {
  reduced_exons <- reduce(x, ignore.strand = TRUE)
  # Add back the gene_id as a metadata column
  mcols(reduced_exons)$gene_id <- mcols(x)$gene_id[1]
  return(reduced_exons)
})

# Combine the list of collapsed exons back into a single GRanges object
all(sapply(collapsed_exons_list, is, "GRanges"))
collapsed_exons <- unlist(GRangesList(collapsed_exons_list))


# Make sure gene_id is a metadata column, not in the row names
names(collapsed_exons) <- NULL  # Remove names (which contain gene_id)

# Move gene_id to the metadata (mcols) as a column
collapsed_exons$gene_id <- mcols(collapsed_exons)$gene_id

# Reorder columns for display (optional)
mcols(collapsed_exons) <- mcols(collapsed_exons)[, "gene_id"]

# View the result
collapsed_exons
saveRDS(collapsed_exons, "collapsed_exons_to_work.RDS")



# Convert your SNP data to a GRanges object
snp_gr <- GRanges(seqnames = snp_data$V1, #join_overlap_inner object
                  ranges = IRanges(start = snp_data$V2, end = snp_data$V2 +1))


#snp_gr <- GRanges(seqnames = snp_data$V1, #subsetOverlap approach
#                  ranges = IRanges(start = snp_data$V2, end = snp_data$V2 + 1,
#                                   names = paste0("rsid:", seq_len(nrow(snp_data)))))


#metadata column                                   
snp_gr$snp_pos <- start(snp_gr)

# Set a 500 kb window around each SNP
#window_size <- 500000
start(snp_gr) <- start(snp_gr) - 500000
end(snp_gr) <- end(snp_gr) + 500000


# =====================================subsetOverlap approach ============================================
#https://support.bioconductor.org/p/54470/
names(snp_gr)
ranges <- subsetByOverlaps(exon_data, snp_gr)
hits <- findOverlaps(exon_data, snp_gr) #90314 - same as join_inner_overlap
rsid <- CharacterList(split(names(snp_gr)[subjectHits(hits)],queryHits(hits)))
snp_pos <- CharacterList(split(snp_gr$snp_pos[subjectHits(hits)],queryHits(hits)))
mcols(ranges) <- DataFrame(mcols(ranges), rsid, snp_pos)



# Convert ranges to data.frame, unnest the CharacterList column, and group by gene_id
target_regions_sub <- as.data.frame(ranges) %>%
  # Convert CharacterList to a standard character vector
  mutate(snp_pos = sapply(snp_pos, paste, collapse = ",")) %>% #since the snp_pos is CharList - mutate to vector
  group_by(gene_id) %>%
  summarise(start_coords = paste(start, collapse = ";"),
            end_coords = paste(end, collapse = ";"),
            chr = unique(seqnames),
            snp_pos = paste(unique(snp_pos), collapse = ";"))

target_regions_sub_unique = unique(target_regions_sub) #1,371


# Split and expand rows function fro snp_pos column
split_and_expand_snp_pos <- function(df) {
  df %>%
    # Separate the target_start and target_end columns into lists
    mutate(
      snp_list = strsplit(as.character(snp_pos), ",")
    ) %>%
    # Unnest the lists into separate rows
    unnest(cols = snp_list) %>%
    # Convert the columns to numeric
    mutate(
      snp_pos = as.numeric(snp_list)
    ) %>%
    # Select the relevant columns
    select(-snp_list)
}

target_regions_sub_unique_split <- split_and_expand_snp_pos(target_regions_sub_unique)
target_regions_sub_unique_split_nona <- target_regions_sub_unique_split %>% drop_na()
target_regions_sub_unique_split_nona = unique(target_regions_sub_unique_split_nona)


# Merge back with the SNP data
final_data_sub <- merge(target_regions_sub_unique_split_nona,snp_data, by.x = "snp_pos", by.y = "V2")
dim(final_data_sub) #2814
dim(final_data) #2820


final_data_unique_sub = unique(final_data_sub)
dim(final_data_unique_sub) #2331
dim(final_data_unique) #2337

colnames(final_data_unique_sub)


#chr, sno_pos, end, ref, alt, strand, test_name, start-target, end-target
final_data_unique_sub_reorder <- final_data_unique_sub %>%
  mutate(chr = chr,
         start = snp_pos,
         end = snp_pos,
         ref_allele = V4, # Replace with actual reference allele
         alt_allele = V5, # Replace with actual alternate allele
         strand = "+",
         test_name = V7,
         target_start = start_coords,
         target_end = end_coords) %>%
  select(chr,start,end,ref_allele, alt_allele,strand, test_name, target_start, target_end)


# Split and expand rows function
split_and_expand <- function(df) {
  df %>%
    # Separate the target_start and target_end columns into lists
    mutate(
      target_start_list = strsplit(as.character(target_start), ";"),
      target_end_list = strsplit(as.character(target_end), ";")
    ) %>%
    # Unnest the lists into separate rows
    unnest(cols = c(target_start_list, target_end_list)) %>%
    # Convert the columns to numeric
    mutate(
      target_start = as.numeric(target_start_list),
      target_end = as.numeric(target_end_list)
    ) %>%
    # Select the relevant columns
    select(-target_start_list, -target_end_list)
}

# Apply the function
result_df_sub <- split_and_expand(final_data_unique_sub_reorder)

#get unique values
result_df_sub_unique = unique(result_df_sub) #21259 rows

#save file
write.table(result_df_sub_unique, "peaks_rna_subsetOverlap.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#================================= Plyranges =========================================================

#2 approach get overlap - join_overlap_inner
#overlapping_exons_ply <- plyranges::join_overlap_inner(exon_data, snp_gr) 
overlapping_exons_ply_1 <- plyranges::join_overlap_inner(collapsed_exons, snp_gr) #18222


# Group exon coordinates by gene for each SNP - ranges approach
target_regions <- overlapping_exons_ply_1 %>%
  as.data.frame() %>%
  group_by(gene_id) %>%
  summarise(start_coords = paste(start, collapse = ";"),
            end_coords = paste(end, collapse = ";"),chr = seqnames, snp_pos = snp_pos)
            
target_regions_unique = unique(target_regions) #A tibble: 1,705 × 5 #Groups:   gene_id [1,371]
exon_data[exon_data$gene_id == "ENSG00000264357",] #check chr info here versus chr info target regions

# Merge back with the SNP data
final_data <- merge(target_regions_unique,snp_data, by.x = "snp_pos", by.y = "V2")
dim(final_data)

final_data_unique = unique(final_data)
dim(final_data_unique)

colnames(final_data_unique)


#chr, sno_pos, end, ref, alt, strand, test_name, start-target, end-target
final_data_unique_reorder <- final_data_unique %>%
  mutate(chr = chr,
         start = snp_pos,
         end = snp_pos,
         ref_allele = V4, # Replace with actual reference allele
         alt_allele = V5, # Replace with actual alternate allele
         strand = "+",
         test_name = V7,
         target_start = start_coords,
         target_end = end_coords) %>%
  select(chr,start,end,ref_allele, alt_allele,strand, test_name, target_start, target_end)


View(final_data_unique_reorder) #ok - 9columns
write.table(final_data_unique_reorder, "peaks_rna_not_splited_fixed_1.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#======================== additional code for aggregating the data - NOT INTERESTED ANYMORE =================
df_aggregated <- final_data_unique_reorder %>%
  group_by(chr,start,end,ref_allele, alt_allele,strand, test_name) %>%
  summarise(
    target_start = list(unlist(strsplit(target_start, ";"))),
    target_end = list(unlist(strsplit(target_end, ";"))),
    .groups = 'drop'
  ) %>%
  # Check if lengths match before combining
  mutate(
    peak_start_len = lengths(target_start),
    peak_end_len = lengths(target_end)
  ) %>%
  # Keep only rows where lengths match
  filter(peak_start_len == peak_end_len) %>%
  # Combine back into a single string for both peak_start and peak_end
  mutate(
    target_start = sapply(target_start, function(x) paste(unique(x), collapse = ";")),
    target_end = sapply(target_end, function(x) paste(unique(x), collapse = ";"))
  )


View(df_aggregated)


# =============Aggregate peak_start and peak_end for the same snp - NOT INTERESTED ANYMORE ============

#1 - checl lengh coord
df <- final_data_unique_reorder
df$target_start_count <- sapply(strsplit(df$target_start, ";"), length)
df$target_end_count <- sapply(strsplit(df$target_end, ";"), length)

View(df)


df_combined_1 <- df %>%
  group_by(chr, start, end, ref_allele, alt_allele, strand, test_name) %>%
  summarise(
    # Split and pair up target_start and target_end together
    paired_values = list(unique(data.frame(
      target_start = unlist(strsplit(target_start, ";")),
      target_end = unlist(strsplit(target_end, ";"))
    ))),
    .groups = 'drop'
  ) %>%
  # Recombine paired values back into strings
  mutate(
    target_start = sapply(paired_values, function(x) paste(x$target_start, collapse = ";")),
    target_end = sapply(paired_values, function(x) paste(x$target_end, collapse = ";"))
  ) %>%
  select(-paired_values) 

df_combined_1$target_start_count <- sapply(strsplit(df_combined$target_start, ";"), length)
df_combined_1$target_end_count <- sapply(strsplit(df_combined$target_end, ";"), length)

View(df_combined_1)



# Ensure that peak_start and peak_end are correctly paired
df_aggregated <- final_data_unique_reorder %>%
  group_by(chr,start,end,ref_allele, alt_allele,strand, test_name) %>%
  summarise(
    target_start = list(unlist(strsplit(target_start, ";"))),
    target_end = list(unlist(strsplit(target_end, ";"))),
    .groups = 'drop'
  ) %>%
  # Check if lengths match before combining
  mutate(
    peak_start_len = lengths(target_start),
    peak_end_len = lengths(target_end)
  ) %>%
  # Keep only rows where lengths match
  filter(peak_start_len == peak_end_len) %>%
  # Combine back into a single string for both peak_start and peak_end
  mutate(
    target_start = sapply(target_start, function(x) paste(unique(x), collapse = ";")),
    target_end = sapply(target_end, function(x) paste(unique(x), collapse = ";"))
  )

# Check for mismatches (Optional debugging step)
mismatches <- final_data_unique_reorder %>%
  group_by(chr,start,end,ref_allele, alt_allele,strand, test_name) %>%
  summarise(
    peak_start = list(unlist(strsplit(target_start, ";"))),
    peak_end = list(unlist(strsplit(target_end, ";"))),
    .groups = 'drop'
  ) %>%
  mutate(
    peak_start_len = lengths(peak_start),
    peak_end_len = lengths(peak_end)
  ) %>%
  filter(peak_start_len != peak_end_len)

# View the aggregated dataframe
df_aggregated_to_save <- df_aggregated[, -c(ncol(df_aggregated)-1, ncol(df_aggregated))]


write.table(df_aggregated_to_save, "peaks_rna_not_splited_fixed.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##splited approach - join overlap
# Apply the function - Split and expand rows
result_df <- split_and_expand(final_data_unique_reorder)

#Get unique rows
result_df_unique = unique(result_df)


# Write the file in the required format
write.table(result_df_unique, "peaks_rna.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



# =======================  Checking rows not included in results dataframes ========================================

# Find rows in result_df_unique that are not in result_df_sub_unique - 153
not_in_sub_unique <- anti_join(result_df_unique, result_df_sub_unique, 
                               by = c("chr", "start", "end", "ref_allele", "alt_allele", "strand", "test_name", "target_start", "target_end"))

write.table(not_in_sub_unique, "peaks_rna_df_unique_not_sub_unique.txt", sep = "\t", quote = FALSE, row.names = FALSE)



# Find rows in result_df_sub_unique that are not in result_df_unique - 0
not_in_unique <- anti_join(result_df_sub_unique, result_df_unique, 
                           by = c("chr", "start", "end", "ref_allele", "alt_allele", "strand", "test_name", "target_start", "target_end"))
