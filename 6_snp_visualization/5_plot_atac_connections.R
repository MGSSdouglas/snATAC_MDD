#Generalize the function above to any two GRanges regions provided to be highlighted #1000
plot_ATAC_connections_two_highlighted_regions <- function(proj=proj, region_1, region_2, left_extension=NULL, 
                                                            right_extension = NULL, highlight_extend=NULL, corCutOff=NULL, chr,
                                                            groupBy="ClustersMapped", useGroups=NULL,
                                                            baseSize=10, facetbaseSize=10, prefix,col=NULL,
                                                            peak_to_gene_loops, #pass this in as an argument
                                                            external_features, external_feature_names, #add in external_feature_palettes
                                                            external_loops, external_loop_names, celltype=NULL) { #add in external_loop_palettes eventually

    #INPUTS: - proj: the name of the ArchR project
    #        - region_1: a numerical position corresponding to the first region to highlight
    #        - region_2: a numerical position corresponding to the second region to highlight
    #        - extension: the width (in base pairs) of the region surrounding the LD block/eQTL that 
    #                     we want to specify (user-specified)
    #        - highlight_extend: the width of the highlight that we want to use to visualize the eQTL/source SNP
    #        - corCutOff: the correlation cutoff for generating peak-to-peak/peak-to-gene
    #        - chr: the chromosome that both loci are on
    #        - groupBy: the feature name of the ArchR project that will be used to generate the bulkTracks
    #        - useGroups: the names of the groups from groupBy that we want to display in bulkTracks
    #        - baseSize/facetbaseSize: font sizes 
    #        - prefix: a string used to identify the plot that has been generated
    #        - celltype: a character used to identify cell-type specific peaks to plot

    #Create our highlight coordinates
    highlight_1 <- region_1
    start(highlight_1) <- start(region_1) - highlight_extend
    end(highlight_1) <- end(region_1) + highlight_extend

    highlight_2 <- region_2
    start(highlight_2) <- start(region_2) - highlight_extend
    end(highlight_2) <- end(region_2) + highlight_extend

    #Make a list of the upper and lower bounds of the positions to include in the plot (before extension)
    positions_to_include <- c(start(highlight_1), end(highlight_1), start(highlight_2), end(highlight_2))
    
    #compute outside of this function as a global variable 

    region <- GRanges(data.frame(chr=chr, start=min(positions_to_include)-left_extension, end=max(positions_to_include)+right_extension))

    #subset p2g for plotting max gene correlation 
    #peak_to_gene_loops= peak_to_gene_loops[peak_to_gene_loops$peakName %in% region_1 & peak_to_gene_loops$Correlation==max(Correlation),]
    #cell-type and cluster peaks
    gc_sub <- readRDS("~/projects/rrg-gturecki/snp_visualization/data/IterativelyRemovedPeaks_sub.rds")
    gc_broad <- readRDS("~/projects/rrg-gturecki/snp_visualization/data/IterativelyRemovedPeaks.rds")     
    #get marker and DAR peak for the cell-type
    marker_all <- read.csv("~/projects/rrg-gturecki/snp_visualization/data/broad_sub_marker_peaks.csv")
    dar_all <- read.csv("~/projects/rrg-gturecki/snp_visualization/data/broad_sub_dar.csv")
    #filter by celltype
    marker <- marker_all[marker_all$celltype %in% celltype,]
    dar <- dar_all[dar_all$cluster_id %in% celltype,]
    #filter to relevant peaks
    #allpeaks <- getPeakSet(proj)
    allpeaks <- c(gc_sub,gc_broad)
    allpeaks$peaks <- paste(seqnames(allpeaks),ranges(allpeaks),sep="-")
    marker_peaks <- allpeaks[allpeaks$peaks %in% c(marker$gene)]
    dar_peaks <- allpeaks[allpeaks$peaks %in% c(dar$gene)]
    peaks <- c(marker_peaks,dar_peaks)
    print(paste0("No. of candidate peaks: ", length(peaks)))
    names(marker_peaks) <- rep("Markers", length(marker_peaks))
    names(dar_peaks) <- rep("DARs", length(dar_peaks))
    ##Convert data into a GRangesList for visualization
    peak_list <- GRangesList(marker_peaks, dar_peaks)
    names(peak_list) <- c("Markers","DARs")
     
#pal=c("#145F04", "#E3D2A7", "#b36609" , "#4853A2" , "#1792AA" , "#C60E9D" ,"#f08080")
#pal=c("salmon4","tan")

    p1 <- plotBrowserTrackWithExternalData( #new function in ArchRBrowserWithExternalData.R
        ArchRProj = proj, 
        baseSize=baseSize,
        facetbaseSize=facetbaseSize,
        features=peak_list,   #peak_list        #getPeakSet(ArchRProj=proj) #can enter the resized GRanges object instead 
        external_features=external_features, #plot LD/eQTL SNPs on different tracks
        external_feature_names=external_feature_names,
        loops = peak_to_gene_loops, #plot P2G by default - can this be changed?
        external_loops=external_loops, #doesn't need to be converted to a list object
        external_loop_names=external_loop_names,
        #groupBy = "predictedBroad_Un",  #change from Clusters
        groupBy = groupBy,pal=col,
        useGroups = useGroups,
        region = region,
        highlight = c(highlight_1, highlight_2), 
        logFile = createLogFile("plotBrowserTrack", logDir="~/scratch/SNP_plots") 
    )

    #plot_list <- c(plot_list, p1)

    #grid::grid.newpage()
    #grid::grid.draw(p1)

    region_1_string <- paste(as.data.frame(seqnames(region_1))[1,1], start(region_1), end(region_1), sep="_")
    region_2_string <- paste(as.data.frame(seqnames(region_2))[1,1], start(region_2), end(region_2), sep="_")
    ggsave(paste0("~/scratch/final_subcluster_plot/", prefix, "_",region_1_string, "_", region_2_string,
                        "_", corCutOff, ".pdf"), p1)

}
