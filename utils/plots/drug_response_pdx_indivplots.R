#' mRECIST labeling
#' 
mrecist_labs <- c(
    'CR' = 'Complete\nResponse',
    'PR' = 'Partial\nResponse',
    'SD' = 'Stable\nDisease',
    'PD' = 'Progressive\nDisease'
)

#' Plot waterfall plot of ARCHE score coloured by mRECIST 
#'
#' Creates waterfall plot of PDX ARCHE scores coloured by mRECIST
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of waterfall plot.
#' 
waterfall_mRECIST <- function(df, arche, drug, plot.indiv = FALSE) {

    # create ranking order
    df <- df[order(df[[arche]], decreasing = T),]
    df$rank <- 1:nrow(df)

    # obtain min and max values
    y <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create waterfall plot coloured by mRESCIST
    p <- ggplot(df, aes(x = rank, y = .data[[arche]], fill = mRECIST)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_manual(values = mrecist_pal, labels = mrecist_labs) +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.8, 'cm'), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            legend.position = "none"
        ) +
        ylim(c(-y, y)) + labs(x = "PDX Model", y = paste(arche, "Score"), fill = "mRECIST") 
    
    # print out waterfall plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/mRECIST_waterfall_", arche, "_", drug, ".png")
        png(path, width=175, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}


#' Plot average ARCHE score for each mRECIST 
#'
#' Creates bar plot of average ARCHE score for each mRECIST
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of bar plot.
#' 
avg_ARCHE_mRECIST <- function(df, arche, drug, plot.indiv = F) {

    # compute average ARCHE score for each mRECIST category
    df <- df %>%
        group_by(mRECIST) %>%
        summarize(
            mean_score = mean(.data[[arche]], na.rm = TRUE),
            sd_score = sd(.data[[arche]], na.rm = TRUE),
            .groups = "drop"
        )
    df$mRECIST <- factor(df$mRECIST, levels = names(mrecist_pal))

    # obtain min and max values
    y <- max(abs(min(df$mean_score)), max(df$mean_score)) |> ceiling()
    
    # create plot
    p <- ggplot(df, aes(x = mRECIST, y = mean_score)) +
        geom_bar(stat = "identity", color = "black", aes(fill = mRECIST)) + ylim(c(-y, y)) +
        geom_hline(yintercept = 0) +
        geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), width = 0.2) +
        scale_fill_manual(values = mrecist_pal) +
        scale_x_discrete(labels = mrecist_labs) +
        labs(x = "\nmRECIST", y = paste(arche, "Score")) +
        theme_classic() + 
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.position = "none")

    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/mRECIST_avg_", arche, "_", drug, ".png")
        png(path, width = 4, height = 5, res = 600, units = "in")
        print(p)
        dev.off()
    }
    return(p)

}

#' Plot waterfall plot of ARCHE score coloured by ARCHE predictive accuracy (compared to mRECIST outcome) 
#'
#' Creates waterfall plot of PDX ARCHE scores coloured by predictive accuracy
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of waterfall plot.
#' 
waterfall_mRECIST_accuracy <- function(df, arche, drug, plot.indiv = F) {

    # create ranking order
    df <- df[order(df[[arche]], decreasing = T),]
    df$rank <- 1:nrow(df)

    # binarize response
    df$predicted <- ifelse(df[[arche]] >= 0, 1, 0)
    df$outcome <- ifelse(df$mRECIST %in% c('CR', 'PR'), 1, 0)

    # determine incorrect predictions
    df$accuracy <- ifelse(df$predicted == df$outcome, 'Correct', 'Incorrect')

    # obtain min and max values
    y <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create plot
    p <- ggplot(df, aes(x = rank, y = .data[[arche]], fill = accuracy)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_manual(values = binary_pal2) +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.8, 'cm'), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            legend.position = "none"
        ) +
        ylim(c(-y, y)) + labs(x = "PDX Model", y = paste(arche, "Score"), fill = "ARCHE Score\nAccuracy") 
    
    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/mRECIST_waterfall_accuracy", arche, "_", drug, ".png")
        png(path, width=175, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}


#' Plot ROC of ARCHE predicted outcome compared to mRECIST 
#'
#' Creates ROC of ARCHE predictive accuracy compared to mRECIST (binarized)
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ROC plot.
#' 
ROC_mRECIST <- function(df, arche, drug, plot.indiv = F) {

    # binarize ARCHE prediction and mRECIST outcome
    predicted <- ifelse(df[[arche]] >= 0, 1, 0)
    actual <- ifelse(df$mRECIST %in% c('CR', 'PR'), 1, 0)

    if (length(unique(predicted)) > 1 && length(unique(actual)) > 1) {
        # compute ARCHE performance
        pred <- prediction(predicted, actual)
        perf <- performance(pred, "tpr", "fpr")

        # create df to plot performance
        toPlot <- data.frame(
            fpr = perf@x.values[[1]],
            tpr = perf@y.values[[1]]
        )

        # create ROC
        p <- ggplot(toPlot, aes(x = fpr, y = tpr)) +
            geom_line(color = "darkgreen", size = 1.2) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
            labs(
                x = "False Positive Rate",
                y = "True Positive Rate"
            ) +
            theme_classic() +
            theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

        # print out plot if needed
        if (plot.indiv == TRUE) {
            path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/mRECIST_waterfall_accuracy", arche, "_", drug, ".png")
            png(path, width=125, height=125, units='mm', res = 600, pointsize=80)
            print(p)
            dev.off()
        }
        return(p)
    } else {
        return(ggplot())
    }
    
}

#' Plot waterfall plot of TR coloured by ARCHE score
#'
#' Creates waterfall plot of PDX TR response coloured by ARCHE scores 
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of waterfall plot.
#' 
waterfall_TR <- function(df, arche, drug, TR, plot.indiv = F) {

    # get label
    label <- switch(
        TR,
        BR_median = "BR",
        BAR_median = "BAR"
    )

    # order samples
    df[[TR]] <- as.numeric(df[[TR]])
    df <- df[order(df[[TR]], decreasing = T),]
    df$rank <- 1:nrow(df)

    # obtain min and max values
    score <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create plot
    p <- ggplot(df, aes(x = rank, y = .data[[TR]], fill = .data[[arche]])) +
        geom_bar(stat = "identity", color = "black") +
        geom_hline(yintercept = 0) +
        scale_fill_gradientn(colors = c("#BC4749", "#F8F1F8", "#077293"), limits = c(-score, score)) +
        theme_classic() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(1, 1),
            legend.justification.inside = c(1, 1),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.5, 'cm'),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        labs(x = "PDX Model", y = paste(drug, label), fill = paste(arche, "\nScore")) 

    # print out waterfall plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/", TR, "_waterfall_", arche, "_", drug, ".png")
        png(path, width=175, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}

#' Plot scatter plot of TR vs ARCHE score
#'
#' Creates scatter plot of PDX TR response vs ARCHE scores 
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if scatter plot should be saved as an individual file.
#' @return ggplot object of scatter plot.
#' 
scatter_TR <- function(df, arche, drug, TR, plot.indiv = F) {

    # get label
    label <- switch(
        TR,
        BR_median = "BR",
        BAR_median = "BAR"
    )

    # compute pearson's correlation
    pc <- cor.test(df[[arche]], df[[TR]], method = "pearson", alternative = "two.sided")
    cor <- round(pc$estimate, 4)
    pval <- round(pc$p.value, 4)

    # create plot
    p <- ggplot(df, aes(x = .data[[arche]], y = .data[[TR]])) +
        geom_smooth(method='lm', formula= y~x, color = "gray") +
        geom_point(size = 3) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = paste(arche, "Score"), y = paste(drug, label)) +
        ggtitle(paste("PCC:", cor, ", pval:", pval))

    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("data/results/figures/4-DrugResponse/PDX/indiv_plots/", TR, "_scatter", arche, "_", drug, ".png")
        png(path, width=125, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}
