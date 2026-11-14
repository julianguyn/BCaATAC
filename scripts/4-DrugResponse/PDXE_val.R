# validate ARCHE associations in PDXE

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("utils/palettes.R")

###########################################################
# Load data
###########################################################

# load ARCHE predictions
scores <- read.csv("data/results/data/1-Signatures/geneSurrogates/pdxe/rna_ARCHE_predictions.csv")

# load drug response
sen <- read.csv("data/procdata/PDXs/pdxe/sen.csv")
sen$mRECIST <- factor(sen$mRECIST, levels = names(mrecist_pal))


###########################################################
# Format drug response associations
###########################################################

# TODO: subset for drug (if more drugs added)

# add ARCHE scores
scores <- scores[match(sen$model.id, scores$sample),]
for (arche in paste0("ARCHE", 1:6)) {
    sen[[arche]] <- scores[[arche]]
}

###########################################################
# Plot drug response associations
###########################################################

plot_pdxe_val <- function(toPlot, arche, drug) {
    res <- cor.test(toPlot[[arche]], toPlot$best.average.response, method = "pearson")
    corr <- round(res$estimate, 2)
    pval <- round(res$p.value, 2)

    toPlot <- toPlot[order(toPlot[[arche]], decreasing = TRUE),]
    toPlot$model.id <- factor(toPlot$model.id, levels = toPlot$model.id)

    p1 <- ggplot(toPlot, aes(x = model.id, y = .data[[arche]], fill = mRECIST)) +
        geom_col(color = "black") +
        scale_fill_manual(values = mrecist_pal) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none"
        ) +
        labs(y = paste(arche, "Score"), x = "PDX Model")

    p2 <- ggplot(toPlot, aes(x = mRECIST, y = .data[[arche]], fill = mRECIST)) +
        stat_summary(fun = mean, geom = "bar", color = "black") +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values = mrecist_pal) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = paste(arche, "Score"))

    p3 <- ggplot(toPlot, aes(x = .data[[arche]], y = best.average.response)) +
        geom_smooth(method = "lm", color = "gray") +
        geom_point() +
        theme_classic() +
        theme(plot.title = element_text(size = 10, hjust = 0.5)) +
        labs(y = paste(drug, "Response (BAR)"), x = paste(arche, "Score")) +
        ggtitle(paste0("PCC: ", corr, " (pval: ", pval, ")"))

    p <- p1 / p2 / p3
    filename <- paste0("data/results/figures/4-DrugResponse/PDX/pdxe/", drug, "_", arche, ".png")
    ggsave(filename, p, width = 3, height = 7)
}

plot_pdxe_val(sen, "ARCHE4", "Binimetinib")