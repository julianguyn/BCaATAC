#' Create panel of plots based on ARCHE association with mRECIST
#'
#' @param subset_df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name
#' @param drug string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @export Panel of four plots.
#' 
assess_ARCHE_mRECIST <- function(subset_df, arche, drug, plot.indiv = FALSE) {

    # subset for variables of interest
    subset_df <- subset_df[,which(colnames(subset_df) %in% c("model_group", "mRECIST", arche))]
    subset_df$mRECIST <- factor(subset_df$mRECIST, levels = c("CR", "PR", "SD", "PD"))
    subset_df <- subset_df[complete.cases(subset_df),]
    #combinations$N.mre[i] <- nrow(subset_df)

    # get plots
    p1 <- waterfall_mRECIST(subset_df, arche, drug)
    p2 <- avg_ARCHE_mRECIST(subset_df, arche, drug)
    p3 <- waterfall_mRECIST_accuracy(subset_df, arche, drug)
    p4 <- ROC_mRECIST(subset_df, arche, drug)
    p <- arrangeGrob(
        p1, p2, p3, p4, ncol = 4, nrow = 2,
        layout_matrix = rbind(c(1,1,2,2), c(3,3,4,4))
    )

    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/", arche, "_", drug, ".png")
        png(path, width=8, height=6, units='in', res = 600, pointsize=80)
        print(grid::grid.draw(p))
        dev.off()
    }
    return(p)
}

#' Create panel of plots based on ARCHE association with selected TR
#'
#' @param subset_df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name
#' @param drug string. Drug name (should match entry under 'drug' column of df)
#' @param TF string. Name of treatment response to use (TR as continuous value)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @export Panel of two plots.
#' 
assess_ARCHE_TR <- function(subset_df, arche, drug, TR, plot.indiv = FALSE) {

    # subset for variables of interest
    subset_df <- subset_df[,which(colnames(subset_df) %in% c("model_group", TR, arche))]
    subset_df[[TR]] <- as.numeric(subset_df[[TR]])
    subset_df <- subset_df[complete.cases(subset_df),]

    if (nrow(subset_df) > 2) {
        # get plots
        p1 <- waterfall_TR(subset_df, arche, drug, TR)
        p2 <- scatter_TR(subset_df, arche, drug, TR)
        p <- ggarrange(p1, p2, ncol = 1, nrow = 2)

        # print out plot if needed
        if (plot.indiv == TRUE) {
            path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/", TR, "/", arche, "_", drug, ".png")
            png(path, width=5, height=6, units='in', res = 600, pointsize=80)
            print(p)
            dev.off()
        }
    } else {
        p <- NA
    }
    return(p)
}

#' Combine ARCHE associations with various treatment responses
#'
#' Combine outputs of assess_ARCHE_mRECIST() and assess_ARCHE_TR()
#' Possible treatment responses (column names):
#' *** mRECIST
#' *** BR_median
#' *** BAR_median
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param label string. Label of ARCHE scores (PDX20k, PDX50k, or PDXall)
#' @param dir string. Parent directory to save results
#' @export Panel of plots.
#' 
assess_ARCHE_PDX <- function(df, label, dir) {

    # initialize dataframe to store results
    combinations <- matrix(data = NA, nrow = length(unique(df$drug)) * 6, ncol = 8) |> as.data.frame()
    colnames(combinations) <- c(
        "ARCHE", "drug", "pair", "N",
        "PC.BR_median", "pval.BR_median",
        "PC.BAR_median", "pval.BAR_median"
    )
    combinations$ARCHE <- rep(paste0("ARCHE", 1:6), length(unique(df$drug)))
    combinations$drug <- rep(unique(df$drug), each = 6)
    combinations$pair <- paste(combinations$ARCHE, combinations$drug, sep = "_")

    # assess ARCHE associations for all available TR
    print(paste("Total number of combinations:", nrow(combinations)))
    for (i in seq_len(nrow(combinations))) {
    #for (i in c(1:3)) {
    #    print(i)
        arche <- combinations$ARCHE[i]
        drug <- combinations$drug[i]
        print(paste(arche, drug))

        subset_df <- df[df$drug == drug,]
        combinations$N[i] <- nrow(subset_df)

        # mRECIST
        p1 <- assess_ARCHE_mRECIST(subset_df, arche, drug)

        # BR Median
        p2 <- assess_ARCHE_TR(subset_df, arche, drug, "BR_median")

        # BAR Median
        p3 <- assess_ARCHE_TR(subset_df, arche, drug, "BAR_median")

        if (nrow(subset_df) > 2) {
            # compute BR pearson's correlation
            pc <- cor.test(subset_df[[arche]], subset_df$BR_median, method = "pearson", alternative = "two.sided")
            combinations$PC.BR_median[i] <- round(pc$estimate, 4)
            combinations$pval.BR_median[i] <- round(pc$p.value, 4)

            # compute BAR pearson's correlation
            pc <- cor.test(subset_df[[arche]], subset_df$BAR_median, method = "pearson", alternative = "two.sided")
            combinations$PC.BAR_median[i] <- round(pc$estimate, 4)
            combinations$pval.BAR_median[i] <- round(pc$p.value, 4)
        }

        # compile plots into panels
        if (!is.na(p2) && !is.na(p3)) {
            filepath <- paste0("data/results/figures/4-DrugResponse/PDX/", dir, "/", label, "/", arche, "_", drug, ".png")
            png(filepath, width=12, height=5, units='in', res = 600, pointsize=80)
            print(grid::grid.draw(
                arrangeGrob(
                    p1, p2, p3, ncol = 4, nrow = 2,
                    layout_matrix = rbind(c(1,1,2,3), c(1,1,2,3))
                ))) # TODO:: prints NULL here, remove somehow
            dev.off()
        } else {
            print("Not enough observations")
        }
    }

    # write out combinations
    filepath <- paste0("data/results/data/4-DrugResponse/PDX/", dir, "_", label, ".csv")
    write.csv(combinations, file = filepath, quote = FALSE, row.names = FALSE)
    return(combinations)
}
