#read in data
gs_final_report <- read.csv("genome_studio_final_report_good_no_head.txt", sep = "\t")
gs_final_report <- gs_final_report[gs_final_report$Sample.Name != "NA10837", ]

#read in genotyping plate metadata
metadata <- read.csv("sc_genotyping_combined_plates_metadata.csv")
metadata$Gender.. <- stringr::str_trim(metadata$Gender..)
metadata$Experimental.Group.. <- stringr::str_trim(metadata$Experimental.Group..)

#read in sample table from genome studio
gs_samples_table <- read.csv("gs_samples_table.txt" ,sep = "\t")
gs_samples_table$Sample.ID <- unlist(stringi::stri_split_fixed(gs_samples_table$Sample.ID, pattern = "#"))[c(FALSE, TRUE, FALSE, FALSE)]
gs_samples_table <- gs_samples_table[gs_samples_table$Sample.ID != "NA10837", ]

#create lgen file
lgen_file <- gs_final_report[, c("Sample.Index", "Sample.Name", "SNP.Name", "Allele1...Plus", "Allele2...Plus")]
lgen_file$Allele1...Plus[lgen_file$Allele1...Plus == "-"] <- 0
lgen_file$Allele2...Plus[lgen_file$Allele2...Plus == "-"] <- 0

#create map file
map_file <- data.frame(gs_final_report$Chr, gs_final_report$SNP.Name, rep(0, times = nrow(gs_final_report)), gs_final_report$Position)

#create fam file
fam_file <- data.frame(gs_samples_table$Index[match(metadata$Sample.Name.., gs_samples_table$Sample.ID)],
                         metadata$Sample.Name..,
                         rep(0, nrow(metadata)),
                         rep(0, nrow(metadata)),
                         ifelse(metadata$Gender.. == "Male", 1, 2),
                         ifelse(metadata$Experimental.Group.. == "Control", 1, 2))

#write files
write.table(lgen_file, "plus_minus_report.lgen", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(map_file, "plus_minus_report.map", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(fam_file, "plus_minus_report.fam", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
