get_hiC_fragments <- function(neu_neg_path, neu_pos_path) {

    non_neuronal_hiC_dataframe <- as.data.frame(read_tsv(
                                                neu_neg_path))

    neuronal_hiC_dataframe <- as.data.frame(read_tsv(
                                            neu_pos_path))
    #Convert to GRanges loop format by getting the midpoint of the two interacting fragments
    #and sorting so that the start position is before the end position
    colnames(non_neuronal_hiC_dataframe) <- c("chr", "start_tss", "end_tss",
                                            "start_interaction", "end_interaction")
    colnames(neuronal_hiC_dataframe) <- c("chr", "start_tss", "end_tss",
                                            "start_interaction", "end_interaction")

    #Convert all HiC fragments into a peak track
    all_non_neuronal_fragments <- rbind(as.matrix(non_neuronal_hiC_dataframe[,c("chr", "start_tss", "end_tss")]),
                                    as.matrix(non_neuronal_hiC_dataframe[,c("chr","start_interaction","end_interaction")]))
    #all_non_neuronal_fragments <- as.matrix(non_neuronal_hiC_dataframe[,c("chr","start_interaction","end_interaction")])

    all_non_neuronal_fragments <- as.data.frame(all_non_neuronal_fragments)
    #sort so that start is always before end
    all_non_neuronal_fragments <- data.frame(all_non_neuronal_fragments[,1], t(apply(all_non_neuronal_fragments[,c(2:3)], 1, sort)))
    colnames(all_non_neuronal_fragments) <- c("chr", "start", "end")
    #all_non_neuronal_fragments <- GRanges(all_non_neuronal_fragments)
    all_non_neuronal_fragments <-  GRanges(data.frame(seqnames=all_non_neuronal_fragments[,"chr"],start=all_non_neuronal_fragments[,"start"],end=all_non_neuronal_fragments[,"end"]))
    #Repeat for Neuronal
    all_neuronal_fragments <- rbind(as.matrix(neuronal_hiC_dataframe[,c("chr", "start_tss", "end_tss")]),
                                    as.matrix(neuronal_hiC_dataframe[,c("chr","start_interaction","end_interaction")]))
    #all_neuronal_fragments <- as.matrix(neuronal_hiC_dataframe[,c("chr","start_interaction","end_interaction")])

    all_neuronal_fragments <- as.data.frame(all_neuronal_fragments)
    all_neuronal_fragments <- data.frame(all_neuronal_fragments[,1], t(apply(all_neuronal_fragments[,c(2:3)], 1, sort)))
    colnames(all_neuronal_fragments) <- c("chr", "start", "end")
    #all_neuronal_fragments <- GRanges(all_neuronal_fragments)
    all_neuronal_fragments <-  GRanges(data.frame(seqnames=all_neuronal_fragments[,"chr"],start=all_neuronal_fragments[,"start"],end=all_neuronal_fragments[,"end"]))

    #Combine Neuronal/Non-Neuronal Data into Single Tracks and then Colour
    #all_hiC_fragments <- rbind(as.data.frame(all_non_neuronal_fragments), as.data.frame(all_neuronal_fragments))
    #all_hiC_fragments <- GRanges(all_hiC_fragments)
    #names(all_hiC_fragments) <- c(rep("neg", 2*nrow(non_neuronal_hiC_dataframe)), rep("pos", 2*nrow(neuronal_hiC_dataframe)))
    #return(all_hiC_fragments)

    #this will separate NeuN+/- frags in two tracks
    names(all_non_neuronal_fragments) <- rep("NeuN- HiC",length(all_non_neuronal_fragments))
    names(all_neuronal_fragments) <- rep("NeuN+ HiC",length(all_neuronal_fragments))
    all_hiC_fragments <- GRangesList(all_non_neuronal_fragments,all_neuronal_fragments)
    names(all_hiC_fragments) <- c("NeuN- HiC","NeuN+ HiC")
    return(all_hiC_fragments)

}

get_hiC_loops <- function(neu_neg_path, neu_pos_path) {

    non_neuronal_hiC_dataframe <- as.data.frame(read_tsv(
                                                neu_neg_path))

    neuronal_hiC_dataframe <- as.data.frame(read_tsv(
                                            neu_pos_path))

    #Convert to GRanges loop format by getting the midpoint of the two interacting fragments
    #and sorting so that the start position is before the end position
    colnames(non_neuronal_hiC_dataframe) <- c("chr", "start_tss", "end_tss",
                                            "start_interaction", "end_interaction")

    #Convert all HiC fragments into a peak track
    all_non_neuronal_fragments <- rbind(as.matrix(non_neuronal_hiC_dataframe[,c("chr", "start_tss", "end_tss")]),
                                    as.matrix(non_neuronal_hiC_dataframe[,c("chr","start_interaction","end_interaction")]))

    all_non_neuronal_fragments <- as.data.frame(all_non_neuronal_fragments)
    #sort so that start is always before end
    all_non_neuronal_fragments <- data.frame(all_non_neuronal_fragments[,1], t(apply(all_non_neuronal_fragments[,c(2:3)], 1, sort)))
    colnames(all_non_neuronal_fragments) <- c("chr", "start", "end")
    all_non_neuronal_fragments <- GRanges(all_non_neuronal_fragments)

    non_neuronal_condensed_df <- data.frame(non_neuronal_hiC_dataframe[,1],
                                0.5*(non_neuronal_hiC_dataframe[,2] + non_neuronal_hiC_dataframe[,3]),
                                0.5*(non_neuronal_hiC_dataframe[,4] + non_neuronal_hiC_dataframe[,5]) )

    #sort so that start is always before end
    non_neuronal_condensed_df <- data.frame(non_neuronal_condensed_df[,1], t(apply(non_neuronal_condensed_df[,c(2:3)], 1, sort)))
    colnames(non_neuronal_condensed_df) <- c("chr", "start", "end")
    non_neuronal_hiC_loops <- list()
    non_neuronal_hiC_loops[["Non-Neuronal HiC"]] <- GRanges(non_neuronal_condensed_df, value=rep(0))


    #Repeat for Neuronal
    all_neuronal_fragments <- rbind(as.matrix(neuronal_hiC_dataframe[,c("chr", "start_tss", "end_tss")]),
                                    as.matrix(neuronal_hiC_dataframe[,c("chr","start_interaction","end_interaction")]))

    all_neuronal_fragments <- as.data.frame(all_neuronal_fragments)
    all_neuronal_fragments <- data.frame(all_neuronal_fragments[,1], t(apply(all_neuronal_fragments[,c(2:3)], 1, sort)))
    colnames(all_neuronal_fragments) <- c("chr", "start", "end")
    all_neuronal_fragments <- GRanges(all_neuronal_fragments)

    neuronal_condensed_df <- data.frame(neuronal_hiC_dataframe[,1],
                                0.5*(neuronal_hiC_dataframe[,2] + neuronal_hiC_dataframe[,3]),
                                0.5*(neuronal_hiC_dataframe[,4] + neuronal_hiC_dataframe[,5]) )

    neuronal_condensed_df <- data.frame(neuronal_condensed_df[,1], t(apply(neuronal_condensed_df[,c(2:3)], 1, sort)))
    colnames(neuronal_condensed_df) <- c("chr", "start", "end")
    neuronal_hiC_loops <- list()
    neuronal_hiC_loops[["Neuronal HiC"]] <- GRanges(neuronal_condensed_df, value=rep(1))
    neuronal_hiC_loops[["Neuronal HiC"]] 


    #Combine Neuronal/Non-Neuronal Data into Single Tracks and then Colour
    all_hiC_fragments <- rbind(as.data.frame(all_non_neuronal_fragments), as.data.frame(all_neuronal_fragments))
    all_hiC_fragments <- GRanges(all_hiC_fragments)
    names(all_hiC_fragments) <- c(rep("non-neuronal", 2*nrow(non_neuronal_hiC_dataframe)), 
                                rep("neuronal", 2*nrow(neuronal_hiC_dataframe)))

    all_hiC_loops <- rbind(non_neuronal_condensed_df, neuronal_condensed_df)
    all_hiC_loops <- GRanges(all_hiC_loops, 
                            value=c(rep(0, nrow(non_neuronal_condensed_df)), rep(1, nrow(neuronal_condensed_df))))
    names(all_hiC_loops) <- c(rep("neg", nrow(non_neuronal_condensed_df)), rep("pos", nrow(neuronal_condensed_df)))
    hiC_loops <- list()
    hiC_loops[["HiC Loops"]] <- all_hiC_loops
    return(hiC_loops)

}


