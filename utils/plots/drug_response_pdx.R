#' Create panel of plots based on ARCHE association with mRECIST
#'
#' @param subset_df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name
#' @param drug string. Drug name (should match entry under 'drug' column of df)
#' @param combinations dataframe. To store results
#' @param i int. Row in combinations to store results
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @export Panel of four plots.
#' 
assess_ARCHE_mRECIST <- function(subset_df, arche, drug, combinations, i, plot.indiv = FALSE) {

    # subset for variables of interest
    subset_df <- subset_df[,which(colnames(subset_df) %in% c("model_group", "mRECIST", arche))]
    subset_df$mRECIST <- factor(subset_df$mRECIST, levels = c("CR", "PR", "SD", "PD"))
    subset_df <- subset_df[complete.cases(subset_df),]
    combinations$N.mre[i] <- nrow(subset_df)

    # get plots
    p1 <- waterfall_mRECIST(subset_df, arche, drug)
    p2 <- avg_ARCHE_mRECIST(subset_df, arche, drug)
    p3 <- waterfall_mRECIST_accuracy(subset_df, arche, drug, combinations, i)
    p4 <- ROC_mRECIST(subset_df, arche, drug)
    p <- arrangeGrob(p1, p2, p3, p4, ncol = 5, nrow = 2,
        layout_matrix = rbind(c(1,1,1,2,2), 
                              c(3,3,3,4,4)))

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
#' @param combinations dataframe. To store results
#' @param i int. Row in combinations to store results
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @export Panel of two plots.
#' 
assess_ARCHE_TR <- function(subset_df, arche, drug, TR, combinations, i, plot.indiv = FALSE) {

    # subset for variables of interest
    subset_df <- subset_df[,which(colnames(subset_df) %in% c("model_group", TR, arche))]
    subset_df[[TR]] <- as.numeric(subset_df[[TR]])
    subset_df <- subset_df[complete.cases(subset_df),]

    # save number of TR experiments
    if (TR == "AUC") {combinations$N.auc[i] <- nrow(subset_df)}
    if (TR == "slope") {combinations$N.slope[i] <- nrow(subset_df)}
    if (TR == "BR_median") {combinations$N.BR_median[i] <- nrow(subset_df)}
    if (TR == "BAR_median") {combinations$N.BAR_median[i] <- nrow(subset_df)}

    if (nrow(subset_df) > 0) {
        # get plots
        p1 <- waterfall_TR(subset_df, arche, drug, TR)
        p2 <- scatter_TR(subset_df, arche, drug, TR, combinations, i)
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
#' *** AUC
#' *** mRECIST
#' *** slope
#' *** BR_median
#' *** BAR_median
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param xeva string. Label of xeva version
#' @export Panel of plots.
#' 
assess_ARCHE_PDX <- function(df, xeva) {

    # initialize TR variables
    auc <- ifelse("AUC" %in% colnames(df), TRUE, FALSE)
    mre <- ifelse("mRECIST" %in% colnames(df), TRUE, FALSE)
    slp <- ifelse("slope" %in% colnames(df), TRUE, FALSE)
    br <- ifelse("BR_median" %in% colnames(df), TRUE, FALSE)
    bar <- ifelse("BAR_median" %in% colnames(df), TRUE, FALSE)

    # initialize dataframe to store results
    combinations <- matrix(data = NA, nrow = length(unique(df$drug)) * 6, ncol = 18) |> as.data.frame()
    colnames(combinations) <- c(
        "ARCHE", "drug", "pair", "N", 
        "PC.auc", "p.auc", "N.auc",
        "AC.mre", "N.mre",
        "PC.slope", "p.slope", "N.slope",
        "PC.BR_median", "p.BR_median", "N.BR_median",
        "PC.BAR_median", "p.BAR_median", "N.BAR_median"
    )
    combinations$ARCHE <- rep(paste0("ARCHE", 1:6), length(unique(df$drug)))
    combinations$drug <- rep(unique(df$drug), each = 6)
    combinations$pair <- paste(combinations$ARCHE, combinations$drug, sep = "_")

    # assess ARCHE associations for all available TR
    #for (i in seq_len(nrow(combinations))) {
    for (i in c(1:3)) {
        print(i)
        arche <- combinations$ARCHE[i]
        drug <- combinations$drug[i]

        subset_df <- df[df$drug == drug,]
        combinations$N[i] <- nrow(subset_df)

        # AUC
        if (auc == TRUE) {
            print("auc")
            p1 <- assess_ARCHE_TR(subset_df, arche, drug, "AUC", combinations, i)
        }
        # mRECIST
        if (mre == TRUE) {
            print("mre")
            p2 <- assess_ARCHE_mRECIST(subset_df, arche, drug, combinations, i)
        }
        # Slope
        if (slp == TRUE) {
            print("slp")
            p3 <- assess_ARCHE_TR(subset_df, arche, drug, "slope", combinations, i)
            if (is.na(p3)) {slp <- FALSE}
        }
        # BR Median
        if (br == TRUE) {
            #print("br")
            p4 <- assess_ARCHE_TR(subset_df, arche, drug, "BR_median", combinations, i)
            if (is.na(p4)) {br <- FALSE}
        }
        # BAR Median
        if (bar == TRUE) {
            #print("bar")
            p5 <- assess_ARCHE_TR(subset_df, arche, drug, "BAR_median", combinations, i)
            if (is.na(p5)) {bar <- FALSE}
        }

        # compile plots into panels (currently 3 possible combinations)
        if (auc == TRUE & mre == TRUE & slp == FALSE & br == FALSE & bar == FALSE) {
            path = paste0("data/results/figures/4-DrugResponse/PDX/", xeva, "/", arche, "_", drug, ".png")
            png(path, width=10, height=6, units='in', res = 600, pointsize=80)
            print(grid::grid.draw(arrangeGrob(p1, p2, ncol = 3, nrow = 2,
                              layout_matrix = rbind(c(1,2,2), c(1,2,2)))))
            dev.off()
        } else if (auc == TRUE & mre == TRUE & slp == TRUE & br == FALSE & bar == FALSE) {
            path = paste0("data/results/figures/4-DrugResponse/PDX/", xeva, "/", arche, "_", drug, ".png")
            png(path, width=13, height=6, units='in', res = 600, pointsize=80)
            print(grid::grid.draw(arrangeGrob(p1, p2, p3, ncol = 4, nrow = 2,
                              layout_matrix = rbind(c(1,2,2,3), c(1,2,2,3)))))
            dev.off()
        } else if (auc == TRUE & mre == TRUE & slp == FALSE & br == TRUE & bar == TRUE) {
            
        } else {
            message("New combination, need to add!")
        }
    }
    return(combinations)
}
