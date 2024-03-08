# source("step17.annotatr_annots.R")

### Step 0: Load the following modules prior to running this script:
# module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
# module load meme/5.4.1
# module load r/4.2.1

library(data.table)
library(annotatr)
library(GenomicRanges)

# prepared finemap ALL table 
cand_snps_table <- fread("./temp/data/candSNP_locs.tsv", header = T)

##### hg38 locations are used #####

# convert SNP information to GRanges
cand_snps <- GRanges(
    seqnames = cand_snps_table$Chromosome,
    ranges = IRanges(start = cand_snps_table$Start, 
        end = cand_snps_table$End, 
        names = cand_snps_table$Name),
    mcols = cand_snps_table$Name
)

# the following is a more minimal set of annotations 
annots <- c("hg38_genes_promoters", "hg38_genes_cds", "hg38_genes_5UTRs", "hg38_genes_exons", "hg38_genes_introns", "hg38_genes_3UTRs", "hg38_genes_intergenic")

# build the annotations of interest using hg38 
annotations = build_annotations(genome = 'hg38', annotations = annots)

# annotate the finemap snps 
cand_snps_annotated <- annotate_regions(
    cand_snps, 
    annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

# convert granges to data.table
table <- as.data.table(cand_snps_annotated)

# save table as tsv to ./temp/data/annotatR_annotated.candSNP_locs.tsv
fwrite(table, file = "./temp/data/annotatR_annotated.candSNP_locs.tsv", sep = "\t", quote = F, row.names = F)

print("Script finished")
