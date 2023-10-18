library(ArchR)
#load
proj <- loadArchRProject("/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered")
genomeAnnotation = getGenomeAnnotation(proj)
geneAnnotation = getGeneAnnotation(proj)
ArchR:::.requirePackage(genomeAnnotation$genome)
ArchR:::.requirePackage("Biostrings",source="bioc")
BSgenome <- eval(parse(text = genomeAnnotation$genome))
BSgenome <- validBSgenome(BSgenome)
promoterRegion = c(2000, 100)
###
setwd("/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered/PeakCalls/ReplicateCalls")
outDir ="/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_BroadPeaks_Split_Filtered/PeakCalls/ReplicateCalls_Relaxed"
#######
outSummit  <- list.files(path="/home/anjali5/scratch/ArchR_outputs/ArchR_BatchTSS_MajorClusters_MetaAdded_FinePeaks_Split_Filtered/PeakCalls/ReplicateCalls/", pattern="*.rds")
outSummit <- do.call("rbind",strsplit(outSummit, "-"))
coverageMetadata <- ArchR:::.getCoverageMetadata(ArchRProj = proj, groupBy = "SubClusters", minCells = 0)
coverageMetadata <- coverageMetadata[match(outSummit[,1], coverageMetadata$Name),]
summitNamesList <- split(coverageMetadata$Name, coverageMetadata$Group)
#make sure they are ordered correctly
summitNamesList <- summitNamesList[unique(coverageMetadata$Group)]
#outSummitList  <- list("Ast" = c(outSummitList[1],outSummitList[2],outSummitList[3],outSummitList[4]), "End" =c(outSummitList[5],outSummitList[6],outSummitList[7],outSummitList[8]), "ExN" = c(outSummitList[9],outSummitList[10],outSummitList[11],outSummitList[12]),"InN" = c(outSummitList[13],outSummitList[14],outSummitList[15],outSummitList[16]),"Mic" = c(outSummitList[17],outSummitList[18],outSummitList[19],outSummitList[20]),"Oli" = c(outSummitList[21],outSummitList[22],outSummitList[23],outSummitList[24]), "OPC" =c(outSummitList[25],outSummitList[26],outSummitList[27],outSummitList[28]))
outSummitList <- lapply(summitNamesList, function(x) {x[x] <- paste0(x,"-summits.rds")})
#make sure these match 
#unique(names(outSummitList))
#unique(names(summitNamesList))
logFile = createLogFile("addReproduciblePeakSet")
logFile="~/projects/def-cnagy/anjali5/ArchRLogs/ArchR-addReproduciblePeakSet-c6059aa3cf5-Date-2022-03-16_Time-14-07-19.log"
groupPeaks <- ArchR:::.safelapply(seq_along(outSummitList), function(i){
  prefix <- sprintf("Rep. Peaks Group (%s of %s) :", i, length(outSummitList))
  ArchR:::.logDiffTime(sprintf("%s Creating Reproducible Peaks", prefix), tstart, verbose = FALSE, logFile = ".")
  peaks <- ArchR:::.identifyReproduciblePeaks(
    summitFiles =  outSummitList[[i]],
    summitNames =  summitNamesList[[i]],
    reproducibility = "2",
    extendSummits = 250,
    blacklist = genomeAnnotation$blacklist,
    prefix = prefix,
    logFile = logFile
  )
  #print(peaks)
  ArchR:::.logDiffTime(sprintf("%s Annotating and Filtering Peaks", prefix), tstart, verbose = FALSE, logFile = logFile)
  peaks <- sort(sortSeqlevels(peaks))
  peaks <- subsetByOverlaps(peaks, genomeAnnotation$chromSizes, type = "within")
  peaks <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion, logFile = logFile)
  peaks <- peaks[which(mcols(peaks)$N < 0.001)] #Remove N Containing Peaks
  peaks <- peaks[order(peaks$groupScoreQuantile, decreasing = TRUE)]
  #peaks <- head(peaks, groupSummary[names(outSummitList)[i],"maxPeaks"])
  peaks <- head(peaks, 300000)
  mcols(peaks)$N <- NULL #Remove N Column
  print(peaks)
  #print(file.path(outDir, paste0(make.names(names(outSummitList)[i]), "-reproduciblePeaks.gr.rds")))
  #saveRDS(peaks, file.path(outDir, paste0(make.names(names(outSummitList)[i]), "-reproduciblePeaks.gr.rds")))
  return(peaks)
}, threads = 10) %>% GRangesList()

names(groupPeaks) <- names(outSummitList)

#Construct Union Peak Set
ArchR:::..logDiffTime("Creating Union Peak Set!", tstart, verbose = verbose, logFile = logFile)
unionPeaks <- unlist(groupPeaks)
unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)
saveRDS(unionPeaks)
#Summarize Output
peakDF <- lapply(seq_along(groupPeaks), function(x){
  data.frame(Group = names(groupPeaks)[x], table(groupPeaks[[x]]$peakType))
}) %>% Reduce("rbind", .)
peakDF$Group <- paste0(peakDF$Group, "(n = ", tableGroups[peakDF$Group],")")
peakDF <- rbind(data.frame(Group = "UnionPeaks", table(unionPeaks$peakType)), peakDF)
peakDF$Freq <- peakDF$Freq / 1000
metadata(unionPeaks)$PeakCallSummary <- peakDF

#Merge by Reduce and annotate
#combined.peaks <- Signac::reduce(x = c(a,b,c,d,e,f,g))
#combined.peaks <- ArchR:::.fastAnnoPeaks(combined.peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion, logFile = logFile)



