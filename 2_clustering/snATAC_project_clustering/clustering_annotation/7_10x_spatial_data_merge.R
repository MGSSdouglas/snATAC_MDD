library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spatialLIBD)

#http://research.libd.org/spatialLIBD/articles/spatialLIBD.html
#get the SingleCellExperiment object for the Lieber Visium daaset 
#data was downloaded using the Globus endpoint provided by https://github.com/LieberInstitute/HumanPilot

sce <- fetch_data(type = 'sce')
male_sections_list <- c("151507", "151508", "151509", "151510")
female_sections_list <- c("151673", "151674", "151675", "151676")

#add metadata
for(section_list in c("male_sections_list", "female_sections_list")){
  print(section_list)
  for(section in sections_list) {
    assign(paste("meta", section, sep = "_"),colData(sce)[colData(sce)$sample_name == section,])
    temp <- Load10X_Spatial(data.dir = paste0("/home/malosree/scratch/spatial_Lieber_data/", section), filename = paste(section,"", sep ="_"))
    temp@meta.data <- cbind(temp@meta.data, get(paste("meta", section, sep="_")))
    assign(paste("section", section, sep = "_"), temp)
  }
  #normalize spatial data
  for(section in ls(pattern = "^section_")){
    assign(section, SCTransform(get(section), assay = "Spatial", verbose = TRUE))
  }
  #set the identities to the layer annotation 
  for(section in ls(pattern = "^section_")){
    temp <- get(section)
    temp$new_layer <- as.character(temp$layer_guess_reordered_short)
    temp$new_layer[is.na(temp$new_layer)]<- rep("N", sum(is.na(temp$new_layer)))
    Idents(temp) <- temp$new_layer
    assign(section, temp)
  }
  #merge slices into abjects
  if(section_list=="male_sections_list"){
    merged_sections <- merge(x= section_151507, y = c(section_151508, section_151509, section_151510), add.cell.ids = c("S151507", "S151508", "S151509", "S151510"))
    merged_sections %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30) -> merged_sections 
    saveRDS(merged_sections, "Male_merged_sections_sptial_LIBD.Rds")
  }else{
    merged_sections <- merge(x= section_151673, y = c(section_151674, section_151675, section_151676), add.cell.ids = c("S151673", "S151674", "S151675", "S151676"))
    merged_sections %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30) -> merged_sections 
    saveRDS(merged_sections, "Female_merged_sections_sptial_LIBD.Rds")
  }
}










