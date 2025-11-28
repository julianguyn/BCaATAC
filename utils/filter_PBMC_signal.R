#' Remove PBMC signal from ARCHE signature sets
#' 
#' @param sites_to_filter bed. Sites from differential analysis between bloor and BCa
#' @param filter_label string. Name for outdir
#' @param overlap int. Percentage overlap between sites needed for removal, default = 25%
#' 
filter_PBMC_signal <- function(sites_to_filter, filter_label, overlap = 25) {

    sites_to_filter <- makeGRangesFromDataFrame(sites_to_filter)
    dir <- "data/procdata/ARCHEs/GriffinSites"
    ARCHEs <- list.files(dir, full.names = TRUE)

    # initialize df to store number of sites removed
    toPlot <- data.frame(filename = sub(paste0(dir, "/"), "", sub("\\.txt", "", ARCHEs)))
    toPlot$ARCHE <- sub("_.*", "", toPlot$filename)
    toPlot$filetype <- sub(".*_", "", toPlot$filename)
    toPlot$num_rm <- NA

    # initialize list to store removed sites
    sites_removed <- list()
    
    message(paste("ARCHE filtering for", filter_label))

    for (file in ARCHEs) {

        filename <- sub(paste0(dir, "/"), "", sub("\\.txt", "", file))
        message(paste("***Filtering", filename))
        bed <- fread(file)

        # initialize new bed file
        filtered_bed <- bed
        colnames(filtered_bed) <- c("Chrom", "Start", "End", "position")

        # define region of overlap
        slice <- 500 * overlap/100
        bed$Start <- bed$Start + slice
        bed$End <- bed$End - slice
        arche <- makeGRangesFromDataFrame(bed)

        # find overlaps
        hits <- findOverlaps(arche, sites_to_filter)
        hitsDF <- data.frame(hits)

        if (!nrow(hitsDF) == 0) {

            # save overlaping sites (that are removed)
            to_rm <- filtered_bed[hitsDF$queryHits, ]
            to_rm <- paste(to_rm$Chrom, to_rm$Start, to_rm$End, sep = ":")
            sites_removed[filename] <- list(to_rm)

            # keep only sites not in overlap
            filtered_bed <- filtered_bed[-hitsDF$queryHits, ]
            
            num_rm <- nrow(bed)-nrow(filtered_bed)
            toPlot$num_rm[toPlot$filename == filename] <- num_rm
            #message(paste("Removed", num_rm, "sites"))

        } else {
            sites_removed[filename] <- ""
        }

        # write out new file
        if (!dir.exists(paste0("data/procdata/ARCHEs/FilteredGriffinSites/", filter_label))) {
            dir.create(paste0("data/procdata/ARCHEs/FilteredGriffinSites/", filter_label), recursive = TRUE)
        }
        outdir <- paste0("data/procdata/ARCHEs/FilteredGriffinSites/", filter_label, "/", filename, ".txt")
        write.table(filtered_bed, file = outdir, quote = F, sep = "\t", col.names = T, row.names = F)
    }

    # plot number of sites removed
    toPlot$prop <- toPlot$num_rm / (as.integer(sub("k", "", toPlot$filetype)) * 1000) * 100
    toPlot$label <- paste0(toPlot$num_rm, " (", toPlot$prop, "%)")
    x <- log(max(toPlot$num_rm) + 0.5) 

    p <- ggplot(toPlot, aes(y = ARCHE, x = log(num_rm), fill = ARCHE)) +
        geom_col() +
        geom_text(aes(x = 0.1, label = label), hjust = 0) +
        scale_x_continuous(lim = c(0, x), expand = c(0,0)) +
        facet_wrap(.~filetype) +
        scale_fill_manual(values = ARCHE_pal) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            legend.key.size = unit(0.5, 'cm')
        ) + 
        labs(
            x = "log(Number of Sites Removed)", 
            title = paste0(filter_label, " (Overlap of ", overlap, "%)")
        )
    #message("Saving figure")
    outdir <- paste0("data/results/figures/5-cfDNA/FilteredGriffinSites/", filter_label, "_", overlap, ".png")
    png(outdir, width=7, height=3, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()

    return(sites_removed)
}

#' Helper function to ARCHE specific lists
#' 
get_list <- function(arche) {
    label_10k <- paste(arche, "10k", sep = "_")
    label_20k <- paste(arche, "20k", sep = "_")
    label_50k <- paste(arche, "50k", sep = "_")

    arche_list <- list(
        all_10k = unlist(removed_all[label_10k], use.names = FALSE),
        all_20k = unlist(removed_all[label_20k], use.names = FALSE),
        all_50k = unlist(removed_all[label_50k], use.names = FALSE),
        er_10k = unlist(removed_er[label_10k], use.names = FALSE),
        er_20k = unlist(removed_er[label_20k], use.names = FALSE),
        er_50k = unlist(removed_er[label_50k], use.names = FALSE),
        basal_10k = unlist(removed_basal[label_10k], use.names = FALSE),
        basal_20k = unlist(removed_basal[label_20k], use.names = FALSE),
        basal_50k = unlist(removed_basal[label_50k], use.names = FALSE),
        her2_10k = unlist(removed_her2[label_10k], use.names = FALSE),
        her2_20k = unlist(removed_her2[label_20k], use.names = FALSE),
        her2_50k = unlist(removed_her2[label_50k], use.names = FALSE)
    )

    return(arche_list)
}

#' Function to plot UPSET plots
#' 
plot_sites_UPSET <- function(list, label = NA, arche_list = FALSE) {

    w <- 10
    h <- 6
    if (arche_list == TRUE) {
        label <- list
        list <- get_list(list)
        w <- 7
        h <- 4
    }
    m <- make_comb_mat(list)
    if (arche_list == FALSE) m <- m[comb_size(m) > 10]
    filename <- paste0("data/results/figures/5-cfDNA/FilteredGriffinSites/", label, ".png")
    png(filename, width=w, height=h, units='in', res = 600)
    print(UpSet(m, set_order = names(arche_list), comb_order = order(-comb_size(m)),
            bg_col = "gray", 
            right_annotation = upset_right_annotation(m)))
    dev.off()
}