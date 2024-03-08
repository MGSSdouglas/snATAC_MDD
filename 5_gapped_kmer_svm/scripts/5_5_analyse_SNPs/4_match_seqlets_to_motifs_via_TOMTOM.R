# source("step10.match_seqlets_to_TF_motifs_via_tomtom.v2.R")

### Step 0: Load the following modules prior to running this script:
# module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
# module load meme/5.4.1
# module load r/4.2.1

### Step 1: "memes" library requires the path to "bin" folder of MEME suite installation. 
meme_bin_path <- "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/MPI/gcc9/openmpi4/meme/5.4.1/bin/"
options(meme_bin = meme_bin_path)

### Configuration 

###### configure file paths to motif databases 
motif_db_basepath <- "../data/motif_data"

# CIS-BP_2.0
CIS_BP_2__path <- paste(motif_db_basepath, 
    "motif_databases/CIS-BP_2.00/Homo_sapiens.meme", sep="/")


library(memes)
library(magrittr)
library(universalmotif)
library(data.table)

check_meme_install()

###### load motif databases

# CISBP
CISBP_motifs <- read_meme(CIS_BP_2__path)
# swap `name` and `altname` of each S4 vector in the list 
CISBP_motifs %<>% 
    lapply(function(x) {
        temp <- x@name
        x@name <- x@altname
        x@altname <- temp
        x@extrainfo <- "CIS-BP_2.0"
        x
    })

###### configure output filepaths 

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = paste0(temp_dir, "/data")

# output filepaths
output_dir = "./output"
output_data_dir = paste0(output_dir, "/data")
output_plot_dir = paste0(output_dir, "/plots")

filepath_base <- paste(temp_data_dir, "TFBS", "tomtom_matches", sep="/")
dir.create(filepath_base, recursive=TRUE, showWarnings=FALSE)

###### start analysis

# load pipeline data table (v6)
snp_table <- paste0(temp_data_dir, "/merged_data.v6.ct_independent.csv")
snp_table <- fread(snp_table)
snp_table <- as.data.frame(snp_table)

# find unique analysis cell types 
unique_cell_types <- unique(snp_table$analysis_cell_type)

# iterate over unique cell types
for (cell_type in unique_cell_types) {

    # filter motif_df for cell type based on analysis_cell_type column 
    cell_type_snp_table <- snp_table[snp_table$analysis_cell_type == cell_type,]

    # create valid motifs from seqlets
    cell_type_snp_table$motif <- mapply(function(motif, name, cell_type) create_motif(motif, name=name, altname=cell_type),
        cell_type_snp_table$seqlet, cell_type_snp_table$Name, cell_type_snp_table$analysis_cell_type)

    # create a universalmotif table from snp_table
    cell_type_motif_df <- cell_type_snp_table$motif %>% to_df()
    cell_type_motif_df$rsid <- cell_type_snp_table$Name
    cell_type_motif_df$analysis_cell_type <- cell_type_snp_table$analysis_cell_type

    ###### TOMTOM motif matching using CISBP motifs
    print(paste("(", cell_type, ") ", "Running TOMTOM on CISBP motifs", sep=""))
    # create cell type specific output filepath
    CISBP_unnested_filepath <- paste(filepath_base, "/", cell_type, "__CISBP_unnested.tsv", sep="")
    # run TOMTOM tool
    tomtom_enhanced_motif_df <- runTomTom(cell_type_motif_df, database=CISBP_motifs)
    # unnest tomtom results
    unnested_df <- tomtom_enhanced_motif_df %>% drop_best_match() %>% tidyr::unnest(tomtom) %>% as.data.frame()
    # drop `match_motif` and `motif` columns from data.frame
    unnested_df <- unnested_df[, -which(names(unnested_df) == "match_motif")]
    unnested_df <- unnested_df[, -which(names(unnested_df) == "motif")]
    unnested_df <- unnested_df[, -which(names(unnested_df) == "bkg")]
    # save results to a file
    write.table(unnested_df, file = CISBP_unnested_filepath, sep = "\t", quote = FALSE, row.names = FALSE)

}

print("Script finished")
