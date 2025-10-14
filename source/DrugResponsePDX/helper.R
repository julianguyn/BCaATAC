#' Standardize PDX Naming
#'
#' Reference defined mapping to standardize PDX names.
#' @param models string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 

# define mapping
mapping <- c("X108099P1" = "108099", 
             "X16720" = "16720",
             "X48602" = "48602",
             "X53782" = "53782",
             "BPTO95" = "BPTO.95",
             "INS_B014" = "INSB014",
             "INS_B019" = "INSB019",
             "REFS032" = "REF032",
             "REFS034P1" = "REF034", 
             "REFS038" = "REF038",
             "REFB047" = "REF047",
             "REF047_S26_L002" = "REF047",
             "REFS036P1A" = "REF036",
             "REF_S_038" = "REF038"
)

map_pdx <- function(models) {

    for (i in 1:length(models)) {
        model = models[i]
        if (model %in% names(mapping)) {
            models[i] <- unname(mapping[model])
        }
    }
    return(models)
}


#' Extract PDX ARCHE signature scores
#'
#' Append columns to treatment response matrices with PDX ARCHE signature scores.
#' @param df string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 
get_ARCHE <- function(df) {
    for (i in 1:nrow(df)) {
        sample = df$patient.id[i]
        df$ARCHE1[i] <- scores[rownames(scores) == sample,]$Signature1
        df$ARCHE2[i] <- scores[rownames(scores) == sample,]$Signature2
        df$ARCHE3[i] <- scores[rownames(scores) == sample,]$Signature3
        df$ARCHE4[i] <- scores[rownames(scores) == sample,]$Signature4
        df$ARCHE5[i] <- scores[rownames(scores) == sample,]$Signature5
        df$ARCHE6[i] <- scores[rownames(scores) == sample,]$Signature6
    }
    return(df)
}

#' Plot waterfall plot of ARCHE score coloured by mRECIST 
#'
#' Creates waterfall plot of PDX ARCHE scores coloured by mRECIST
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of waterfall plot.
#' 
waterfall_mRECIST <- function(df, arche, drug, plot.indiv = F) {
 
    # create ranking order
    df <- df[order(df[[arche]], decreasing = T),]
    df$rank <- 1:nrow(df)

    # obtain min and max values
    y <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create waterfall plot coloured by mRESCIST
    p <- ggplot(df, aes(x = rank, y = .data[[arche]], fill = mRECIST)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_manual(values = c("#136F63", "#9DCBBA", "#FFADA1", "#B02E0C"),
                        labels = c('CR' = 'Complete\nResponse', 
                                'PD' = 'Progressive\nDisease', 
                                'SD' = 'Stable\nDisease',
                                'PR' = 'Partial\nResponse')) +
        theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylim(c(-y, y)) + labs(x = "PDX Model", y = paste(arche, "Score"), fill = "mRECIST") 
    
    # print out waterfall plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("DrugResponsePDX/results/figures/indiv_plots/mRECIST_waterfall_", arche, "_", drug, ".png")
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

    # obtain min and max values
    y <- max(abs(min(df$mean_score)), max(df$mean_score)) |> ceiling()
    
    # create plot
    p <- ggplot(df, aes(x = mRECIST, y = mean_score)) +
        geom_bar(stat = "identity", color = "black", aes(fill = mRECIST)) + ylim(c(-y, y)) +
        geom_hline(yintercept = 0) +
        geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), width = 0.2) +
        scale_fill_manual(values = c("#136F63", "#9DCBBA", "#FFADA1", "#B02E0C"),
                            labels = c('CR' = 'Complete\nResponse', 
                                    'PD' = 'Progressive\nDisease', 
                                    'SD' = 'Stable\nDisease',
                                    'PR' = 'Partial\nResponse')) +
        scale_x_discrete(labels = c('CR' = 'Complete\nResponse', 
                                    'PD' = 'Progressive\nDisease', 
                                    'SD' = 'Stable\nDisease',
                                    'PR' = 'Partial\nResponse')) +
        labs(x = "\nmRECIST", y = paste(arche, "Score")) +
        theme_classic() + 
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.position = "none")

    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("DrugResponsePDX/results/figures/indiv_plots/mRECIST_avg_", arche, "_", drug, ".png")
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
    df$accuracy <- ifelse(df$predicted == df$outcome, 'correct', 'incorrect')

    # obtain min and max values
    y <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create plot
    p <- ggplot(df, aes(x = rank, y = .data[[arche]], fill = accuracy)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_manual(values = c("#1E4076", "#A94745")) +
        theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylim(c(-y, y)) + labs(x = "PDX Model", y = paste(arche, "Score"), fill = "ARCHE Score\nAccuracy") 
    
    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("DrugResponsePDX/results/figures/indiv_plots/mRECIST_waterfall_accuracy", arche, "_", drug, ".png")
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
            path = paste0("DrugResponsePDX/results/figures/indiv_plots/mRECIST_waterfall_accuracy", arche, "_", drug, ".png")
            png(path, width=125, height=125, units='mm', res = 600, pointsize=80)
            print(p)
            dev.off()
        }
        return(p)
    } else {
        return(ggplot())
    }
    
}


#' Create panel of plots based on ARCHE association with mRECIST
#'
#' Uses pre-defined functions (above) to create various ARCHE plots based on mRECIST 
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @export Panel of four plots.
#' 
assess_ARCHE_mRECIST <- function(df, drug) {

    # subset for variables of interest
    df <- df[df$drug == drug,]
    if (c('slope') %in% colnames(df)) {
        df <- df[,-which(colnames(df) %in% c('slope', 'AUC'))]
    }
    df$mRECIST <- factor(df$mRECIST, levels = c("CR", "PR", "SD", "PD"))
    df <- df[complete.cases(df),]

    # iterate through each ARCHE
    for (arche in paste0('ARCHE', 1:6)) {
        p1 <- waterfall_mRECIST(df, arche, drug)
        p2 <- avg_ARCHE_mRECIST(df, arche, drug)
        p3 <- waterfall_mRECIST_accuracy(df, arche, drug)
        p4 <- ROC_mRECIST(df, arche, drug)

        # create figure panel
        path = paste0("DrugResponsePDX/results/figures/mRECIST/", arche, "_", drug, ".png")
        png(path, width=8, height=6, units='in', res = 600, pointsize=80)
        print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2))
        dev.off()
    }
}

#' Plot waterfall plot of AUC coloured by ARCHE score
#'
#' Creates waterfall plot of PDX AUC response coloured by ARCHE scores 
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of waterfall plot.
#' 
waterfall_AUC <- function(df, arche, drug, plot.indiv = F) {

    # order samples
    df$AUC <- as.numeric(df$AUC)
    df <- df[order(df$AUC, decreasing = T),]
    df$rank <- 1:nrow(df)

    # obtain min and max values
    score <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create plot
    p <- ggplot(df, aes(x = rank, y = AUC, fill = .data[[arche]])) + 
        geom_bar(stat = "identity", color = "black") + 
        geom_hline(yintercept = 0) +
        scale_fill_gradientn(colors = c("#BC4749", "#F8F1F8", "#077293"), limits = c(-score, score)) +
        theme_classic() + 
        theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(x = "PDX Model", y = "AUC", fill = paste(arche, "\nScore")) 

    # print out waterfall plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("DrugResponsePDX/results/figures/indiv_plots/AUC_waterfall_", arche, "_", drug, ".png")
        png(path, width=175, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}

#' Plot scatter plot of AUC vs ARCHE score
#'
#' Creates scatter plot of PDX AUC response vs ARCHE scores 
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if scatter plot should be saved as an individual file.
#' @return ggplot object of scatter plot.
#' 
scatter_AUC <- function(df, arche, drug, plot.indiv = F) {

    # compute pearson's correlation
    pc <- cor.test(df[[arche]], df$AUC, method = "pearson", alternative = "two.sided")

    # create plot
    p <- ggplot(df, aes(x = .data[[arche]], y = AUC)) +
        geom_point(size = 2) +
        geom_smooth(method='lm', formula= y~x) +
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = paste(arche, "Score")) +
        ggtitle(paste("PCC:", round(pc$estimate, 4), ", pval:", round(pc$p.value, 4)))

    # print out plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("DrugResponsePDX/results/figures/indiv_plots/AUC_scatter", arche, "_", drug, ".png")
        png(path, width=125, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}


#' Create panel of plots based on ARCHE association with AUC
#'
#' Uses pre-defined functions (above) to create various ARCHE plots based on AUC 
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @export Panel of four plots.
#' 
assess_ARCHE_AUC <- function(df, drug) {

    # subset for variables of interest
    df <- df[df$drug == drug,]
    if (c('slope') %in% colnames(df)) {
        df <- df[,-which(colnames(df) %in% c('slope', 'mRECIST'))]
    }
    df$AUC <- as.numeric(df$AUC)
    df <- df[complete.cases(df),]

    # iterate through each ARCHE
    for (arche in paste0('ARCHE', 1:6)) {
        p1 <- waterfall_AUC(df, arche, drug)
        p2 <- scatter_AUC(df, arche, drug)

        # create figure panel
        path = paste0("DrugResponsePDX/results/figures/AUC/", arche, "_", drug, ".png")
        png(path, width=5, height=6, units='in', res = 600, pointsize=80)
        print(ggarrange(p1, p2, ncol = 1, nrow = 2))
        dev.off()
    }
}