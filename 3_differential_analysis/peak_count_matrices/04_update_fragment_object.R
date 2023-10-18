base="/home/anjali5/projects/def-cnagy/Nanuq_Aug21_snATAC/cellranger_v2_output/"
names=c("A1B1","A2B2","A3B3","A4B4","A5B5","A6B6","A7B7","A8B8","A9B9","A10B10","A11B11","A12B12","A13B13","A14B14","A15B15","A16B16","A17B17","A18B18","A19B19","A20B20","A21B21","A22B22","A23B23","A24B24","A25B25","A26B26","A27B27","A28B28","A29B29","A30B30","A31B31","A32B32","A33B33","A34B34","A35B35","A36B36","A37B37","A38B38","A39B39","A40B40","A41B41","A42B42")
new.list=list()
for(name in names){
  tmp = paste0(base,name,"/fragments.tsv.gz")
  new.list <- append(new.list,tmp)
}
for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.list[[i]]) # update path
}
