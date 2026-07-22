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
    subset_df <- subset_df[,which(colnames(subset_df) %in% c("model_group", "mRECIST", "Subtype", arche))]
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
        path <- paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/", arche, "_", drug, ".png")
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
    subset_df <- subset_df[,which(colnames(subset_df) %in% c("model_group", "Subtype", TR, arche))]
    subset_df[[TR]] <- as.numeric(subset_df[[TR]])
    subset_df <- subset_df[complete.cases(subset_df),]

    if (nrow(subset_df) > 2) {
        # get plots
        p1 <- waterfall_TR(subset_df, arche, drug, TR)
        p2 <- scatter_TR(subset_df, arche, drug, TR)
        p <- ggarrange(p1, p2, ncol = 1, nrow = 2)

        # print out plot if needed
        if (plot.indiv == TRUE) {
            path <- paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/", arche, "_", drug, ".png")
            png(path, width=3.5, height=5, units='in', res = 600, pointsize=80)
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
assess_ARCHE_PDX <- function(df, label, dir, meta, plot = TRUE, pcc_thres = 0.4) {

    # initialize dataframe to store results
    combinations <- matrix(data = NA, nrow = length(unique(df$drug)) * 6, ncol = 10) |> as.data.frame()
    colnames(combinations) <- c(
        "ARCHE", "drug", "pair", "N",
        "PC.BR_median", "pval.BR_median",
        "PC.BAR_median", "pval.BAR_median",
        "PDX_subset", "ARCHE_label"
    )
    combinations$ARCHE <- rep(paste0("ARCHE", 1:6), length(unique(df$drug)))
    combinations$drug <- rep(unique(df$drug), each = 6)
    combinations$pair <- paste(combinations$ARCHE, combinations$drug, sep = "_")
    combinations$PDX_subset <- dir
    combinations$ARCHE_label <- label

    # assess ARCHE associations for all available TR
    print(paste("Total number of combinations:", nrow(combinations)))
    for (i in seq_len(nrow(combinations))) {
    #for (i in c(1:3)) {
    #    print(i)
        arche <- combinations$ARCHE[i]
        drug <- combinations$drug[i]
        #print(paste(arche, drug))

        subset_df <- df[df$drug == drug,]
        subset_df <- subset_df[
            complete.cases(subset_df[, c("BR_median", "BAR_median")]),
        ]

        combinations$N[i] <- nrow(subset_df)

        if (plot == TRUE) {
            # mRECIST
            p1 <- assess_ARCHE_mRECIST(subset_df, arche, drug)

            # BR Median
            p2 <- assess_ARCHE_TR(subset_df, arche, drug, "BR_median")

            # BAR Median
            p3 <- assess_ARCHE_TR(subset_df, arche, drug, "BAR_median")
        }

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

        if (plot == TRUE) {
            # compile plots into panels
            if (!is.na(p2) && !is.na(p3) && (abs(combinations$PC.BAR_median[i]) > pcc_thres || abs(combinations$PC.BR_median[i]) > pcc_thres)) {
                filepath <- paste0("data/results/figures/4-DrugResponse/PDX/", dir, "/", label, "/", arche, "_", drug, ".png")
                png(filepath, width=12, height=5, units='in', res = 600, pointsize=80)
                grid::grid.draw(
                    arrangeGrob(
                        p1, p2, p3, ncol = 4, nrow = 2,
                        layout_matrix = rbind(c(1,1,2,3), c(1,1,2,3))
                    )) # TODO:: prints NULL here, remove somehow
                dev.off()
            }
        }
    }

    #if (plot == TRUE) {
    #    # write out combinations
    #    filepath <- paste0("data/results/data/4-DrugResponse/PDX/", dir, "/", label, ".csv")
    #    write.csv(combinations, file = filepath, quote = FALSE, row.names = FALSE)
    #}
    #return(combinations)
}

#' Function to create bubble plots of PDX drug response associations
#' 
#' @param compile data.frame. Compiled from assess_ARCHE_PDX(), rbinded
#' 

plot_compiled_PDX <- function(compile) {

    sig <- compile[which(abs(compile$PC.BAR_median) > 0.4 & compile$pval.BAR_median < 0.1),]
    sig_pairs <- sig$pair
    toPlot <- compile[compile$pair %in% sig_pairs,]
    toPlot$sig <- ifelse(toPlot$pval.BAR_median < 0.1, 'pval < 0.1', 'pval >= 0.1')

    toPlot$Label <- factor(toPlot$Label, levels = c(
        "zscore_T", "zscore_K", "normzscr_T", "normzscr_K",
        "zscore_sumdev_T", "zscore_sumdev_K", "normzscr_sumdev_T", "normzscr_sumdev_K"
        )
    )

    for (arche in paste0("ARCHE", 1:6)) {

        subset <- toPlot[toPlot$ARCHE == arche,]
        subset$score <- subset$Label

        p1 <- ggplot(subset, aes(x = drug, y = score, fill = N)) +
            geom_tile() +
            geom_text(aes(label = N), size = 2.5) +
            scale_fill_gradient(high = "#B6B8D6", low = "#BBDBD1") +
            theme_void() +
            theme(
                axis.text.y = element_text(size = 8, hjust = 1),
                legend.position = "none") +
            ggtitle(arche)

        p2 <- ggplot(subset, aes(x = drug, y = Label, fill = PC.BAR_median, size = -log(pval.BAR_median), shape = sig)) +
            geom_point() +
            geom_text(data = subset(subset, sig == 'pval < 0.1'), aes(label = round(PC.BAR_median, 2)), color = "black", size = 2.5) +
            scale_shape_manual(values = c(21, 24)) +
            scale_size(range = c(2, 12)) +
            scale_fill_gradient2(
                low = "#BC4749",
                high = "#689CB0",
                mid = "#C2BBC9",
                limits = c(-1, 1)
            ) +
            theme_bw() +
            theme(
                legend.key.size = unit(0.3, 'cm'),
                axis.text.x = element_text(size=6, angle=25, hjust=1, vjust=1, margin = margin(t = 3))
            ) 
        
        p <- p1 / p2 + plot_layout(heights = c(3, 6))
        filename <- paste0("data/results/figures/4-DrugResponse/PDX/compiled/", arche, ".png")
        ggsave(filename, p, w = 10, h = 5)
    }

}

#' Plot ARCHE associations for ADCs 
#' 
#' @param scores matrix
#' @param label string folder name
#' @param drug string drug name
#' @param target string Ensembl ID of gene of protein target
#' @param threshold string percentage of quantile to use as lower bound
#' 
plot_ADC_ARCHE <- function(scores, label, drug, target, threshold = "25%") {
    
    # subset rna by gene expression
    subset_rna <- rna[target,]
    threshold <- quantile(as.numeric(subset_rna))[threshold]
    high_exp <- colnames(subset_rna)[subset_rna > threshold]
    low_exp <- colnames(subset_rna)[subset_rna <= threshold]

    # subset ARCHE expression
    toPlot <- scores[scores$drug == drug,]
    high_exp <- toPlot[toPlot$PDX_ID %in% high_exp,]
    low_exp <- toPlot[toPlot$PDX_ID %in% low_exp,]

    assess_ARCHE_PDX(high_exp, label, "adcs", plot = TRUE, pcc_thres = 0)
}