
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(ggplot2)
library(muscat)
library(Seurat)
library(limma)

#script pulls E-values of limma-voom from muscat.

resDS.mod <- function(x, y, bind = c("row", "col"),
    frq = FALSE, cpm = FALSE, digits = 3, sep = "_", ...) {
    
    # check validity of input arguments
    muscat:::.check_sce(x, req_group = TRUE)
    #.check_res(x, y)
    bind <- match.arg(bind)
    if (!is.logical(frq)) 
        muscat:::.check_frq(x, frq)
    stopifnot(is.infinite(digits) || is.numeric(digits) &
        digits > 0 & as.integer(digits) == digits)

    ei <- metadata(x)$experiment_info
    kids <- levels(x$cluster_id)
    
    res <- switch(bind,
        row = {
            bind_rows(lapply(y$table, bind_rows))
        },
        col = {
            ct <- ifelse(!is.null(y$args$contrast), "contrast", "coef")
            cs <- names(y$table)
            res <- lapply(cs, function(c) {
                df <- bind_rows(y$table[[c]])
                df <- select(df, -ct)
                i <- !colnames(df) %in% c("gene", "cluster_id")
                colnames(df)[i] <- paste(colnames(df)[i], c, sep = sep)
                return(df)
            })
            reduce(res, full_join, by = c("gene", "cluster_id"))
        })

    tidy <- function(u, ei, append = "") {
        m1 <- match(ei$sample_id, colnames(u), nomatch = 0)
        m2 <- match(levels(ei$group_id), colnames(u))
        if (all(is.na(m2))) m2 <- 0
        colnames(u)[m1] <- paste0(ei$sample_id, append)[m1 != 0]
        colnames(u)[m2] <- paste0(colnames(u)[m2], append)
        k <- seq_len(ncol(u))[-c(m1, m2)]
        u[, c(k, m1[order(ei$group)], m2)]
    }
    
    # append expression frequencies
    if (is.logical(frq))
        if (frq) frq <- calcExprFreqs(x, ...) else frq <- NULL
    if (!is.null(frq)) {
        frq <- data.frame(
            gene = rep(rownames(x), length(assays(frq))),
            cluster_id = rep(assayNames(frq), each = nrow(x)), 
            do.call("rbind", as.list(assays(frq))),
            row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
        frq <- tidy(frq, ei, append = ".frq")
        res <- inner_join(frq, res, by = c("gene", "cluster_id"))
    }

    # append CPMs
    if (cpm) {
        cpm <- lapply(kids, function(k) {
            if (is.null(y$data[[k]])) 
                return(NULL)
            cpm <- y$data[[k]][[2]] #the line modified
            data.frame(cpm, 
                gene = rownames(cpm),
                cluster_id = k,
                row.names = NULL, 
                check.names = FALSE, 
                stringsAsFactors = FALSE)
        })
        cpm <- bind_rows(cpm)
        cpm <- tidy(cpm, ei, append = ".cpm")
        res <- inner_join(cpm, res, by = c("gene", "cluster_id"))
    }
    mutate_if(res, is.numeric, signif, digits)
}

