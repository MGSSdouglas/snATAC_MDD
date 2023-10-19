#read in plink data
plink_data <- read.table("plink_duplicates_only.traw", header = TRUE)

#make keys to quickly identify duplicated SNPs
keys <- mapply(plink_data$CHR, plink_data$POS, plink_data$ALT, FUN = paste0)

#columns of plink_data that represent samples
columns_to_check <- 7:148

#function to count the total number of matches in SNP call between duplicate SNP
#as well as the total number of comparisons made
check_consistency <- function(data_column, keys) {
  curr_sum <- 0
  total_comparisons <- 0
  for (key in keys) {
		current_snps <- data_column[keys == key]
		current_snps <- current_snps[!is.na(current_snps)]
		if (length(current_snps > 1)) {
		  total_comparisons <- total_comparisons + (length(current_snps) - 1)
		  curr_sum <- curr_sum + sum(duplicated(current_snps))
                  data_column <- data_column[keys != key]
                  keys <- keys[keys != key]
		}
  }
  return(c(curr_sum, total_comparisons))
}

#calculate the percentage of duplicated SNPs that match
match_percents <- vector(mode = "numeric", length = length(columns_to_check))
for(column in columns_to_check) {
  	check_results <- check_consistency(plink_data[ ,column], keys)
  	match_percents[[column - 6]] <- check_results[1] / check_results[2]
}

sample_names <- unlist(lapply(colnames(plink_data[ ,columns_to_check]), FUN = function(x){return(unlist(stringr::str_split(x, "_"))[2])}))

write.csv(data.frame("sample_id" = sample_names, "percent_matches" = match_percents), "duplicate_SNP_consistency_check.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
