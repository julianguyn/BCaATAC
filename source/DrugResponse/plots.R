#' Plot individual associations across PSets
#' 
plot_indivPlot <- function(pair, class) {

    arche <- gsub("_.*", "", pair)
    drug <- gsub(paste0(arche, "_"), "", pair)

    # fix that one drug name
    to_fix = "doxorubicin:navitoclax (2:1 mol/mol)"
    if (drug == to_fix) {
        pair <- "ARCHE1_doxorubicin:navitoclax"
    }

    # set up dataframe to store results
    all_drug_sen <- data.frame(matrix(nrow=0, ncol=3))
    colnames(all_drug_sen) <- c("Sample", "AAC", "PSet")

    # subset to keep only drug of interest
    all_drug_sen <- get_doiAAC(ubr1_sen, drug, all_drug_sen, "UBR1")
    all_drug_sen <- get_doiAAC(ubr2_sen, drug, all_drug_sen, "UBR2")
    all_drug_sen <- get_doiAAC(gray_sen, drug, all_drug_sen, "GRAY")
    all_drug_sen <- get_doiAAC(gcsi_sen, drug, all_drug_sen, "gCSI")
    all_drug_sen <- get_doiAAC(gdsc_sen, drug, all_drug_sen, "GDSC2")
    all_drug_sen <- get_doiAAC(ccle_sen, drug, all_drug_sen, "CCLE")
    all_drug_sen <- get_doiAAC(ctrp_sen, drug, all_drug_sen, "CTRP")

    # get signature scores
    signature_scores <- get_all_scores()
    signature_scores <- signature_scores[rownames(signature_scores) == arche,]
    all_drug_sen$Score <- NA
    for (i in 1:nrow(all_drug_sen)) { 
        all_drug_sen$Score[i] <-  signature_scores[,colnames(signature_scores) == all_drug_sen$Sample[i]]
    }

    # get x axis limits for plotting
    x <- max(max(all_drug_sen$Score), abs(min(all_drug_sen$Score)))

    # set up palette for plotting
    pal <- PSet_pal[names(PSet_pal) %in% unique(all_drug_sen$PSet)]

    png(paste0("DrugResponse/results/figures/",class,"/indivPlots/",pair,".png"), width=160, height=125, units='mm', res = 600, pointsize=80)
    print(ggplot(all_drug_sen, aes(x = Score, y = AAC, fill = PSet)) + 
        geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", se=F, aes(color = PSet)) + 
        scale_fill_manual(values = pal) + 
        scale_color_manual(values = pal) +
        guides(color = 'none', fill = guide_legend(override.aes=list(values = pal[names(pal) %in% unique(all_drug_sen$PSet)], linetype = 0))) +
        theme_classic() + 
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5, size = 16), 
            legend.key.size = unit(0.7, 'cm')
        ) +
        xlim(-x, x) + ylim(0, 1) + 
        labs(
            x = paste0(arche, " Score"), 
            y = paste0(drug, " Response (AAC)"),
            title = paste0(arche, " - ", drug)
        ))
    dev.off()
}

#' Plot Class A associations across PSets
#'
#' @param ARCHE string. ARCHE to subset for
#' @param width int. Width for figure output
#' Colour by concordance of association with ClassA biomarker
#' 
plot_ClassA_allAssociations <- function(toPlot, ARCHE, width) {

    toPlot <- toPlot[toPlot$signature == ARCHE,]
    filename <- paste0("DrugResponse/results/figures/ClassA/supplementary/allAssociations_", ARCHE, ".png")

    png(filename, width = width, height = 5, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = pset, y = pc, fill = ifelse(pc>0, "Positive", "Negative"))) +
        geom_bar(stat="identity", color = "black") +
        geom_text(data = subset(toPlot[toPlot$pc > 0,], FDR < 0.05),
            aes(label = "*", y = pc), 
            vjust = 0, size = 6) +
        geom_text(data = subset(toPlot[toPlot$pc < 0,], FDR < 0.05),
            aes(label = "*", y = pc - 0.17), 
            vjust = 0, size = 6) +
        facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
        scale_fill_manual(values = binary_pal, breaks = c("Positive", "Negative"),
            labels = c("Positive\nAssociation",
                    "Negative\nAssociation")) +
        geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted") +
        geom_hline(yintercept = 0) +
        ylim(c(-1, 1)) + 
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.spacing = unit(0, "lines")
        ) +
        labs(y = "Pearson's Correlation Coefficient", fill = "Direction of\nAssociation", x = "Dataset"))
    dev.off()
}

#' Plot Class A associations as heatmap
#'
#' @param type string. "Multi" for drugs in >1 PSet, "Single" for drugs in only 1 PSet
#' Annotated by ARCHE and BCa relevant drug
#' 
plot_ClassA_heatmap <- function(toPlot, type) {

    # format plot
    toPlot <- toPlot[order(toPlot$pairs),]
    toPlot$pset <- factor(toPlot$pset, levels = names(PSet_pal))
    toPlot$pairs <- factor(toPlot$pairs, levels = unique(toPlot$pairs))

    # arche bounds 
    n_pairs <- toPlot %>%
        group_by(signature) %>%
        summarise(n_pairs = n_distinct(pairs)) %>%
        mutate(cumulative = cumsum(n_pairs))
    bounds <- n_pairs$cumulative+0.5

    # main heatmap (pset~pair)
    p1 <- ggplot(toPlot, aes(x = pairs, y = pset, fill = pc)) + 
        geom_tile(color = 'black') +
        geom_text(data = subset(toPlot, FDR < 0.05),
                aes(label = "*"), 
                vjust = 0.75, size = 4) +
        geom_vline(xintercept = bounds[1:5], color = "gray") +
        scale_x_discrete(labels = function(pairs) gsub(".*_", "", pairs)) +
        scale_fill_gradient2("Pearson's\nCorrelation\nCoefficient", 
                            low = "#BC4749", 
                            high = "#689CB0",
                            mid = "#C2BBC9",
                            limits = c(-1, 1)) +
        theme_void() +
        theme(
            axis.text.y = element_text(size=9, hjust=1, vjust=0.5, margin = margin(r = 6)), 
            axis.text.x = element_text(size=9, angle=90, hjust=1, vjust=0.5, margin = margin(t = 3)),
            legend.title = element_text(size=9),
            axis.ticks = element_line(color = "gray", linewidth = 0.3),
            axis.ticks.length = unit(2, "pt")
            )

    # BCa drug annotation
    p2 <- ggplot(toPlot, aes(x = pairs, y = 1, fill = ifelse(drug %in% bca_drugs, "A", "B"))) + 
        geom_tile(color = "gray") +
        theme_void() + 
        geom_vline(xintercept = bounds[1:5], color = "gray") +
        scale_fill_manual(values = c("#3E517A", "white")) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 2))) + 
        labs(y = "BCaDrug")
    
    # ARCHE annotation
    p3 <- ggplot(toPlot, aes(x = pairs, y = 1, fill = signature)) + 
        geom_tile(color = NA) +
        theme_void() + 
        geom_vline(xintercept = bounds[1:5], color = "gray") +
        scale_fill_manual("ARCHE", values = ARCHE_pal) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 5))) + 
        labs(y = "ARCHE ")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")

    # figure widths
    w <- ifelse(type == "Multi", 18, 21)

    filename <- paste0("DrugResponse/results/figures/ClassA/heatmap", type, ".png")
    png(filename, width = w, height = 3.5, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, p3, l3, l1, ncol = 10, nrow = 17,
        layout_matrix = rbind(c(3,3,3,3,3,3,3,3,3,3,3,3,NA,NA),
                              c(2,2,2,2,2,2,2,2,2,2,2,2,NA,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,4,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA)))
    )
    dev.off()
}

#' Plot Class A biomarker associations
#' 
plot_ClassA_biomarkersAssociations <- function(ClassA) {
    png("DrugResponse/results/figures/ClassA/biomarkersAssociations.png", width = 6, height = 25, res = 600, units = "in")
    print(ggplot(ClassA, aes(x = pc, y = rank)) +
        geom_col(aes(fill = signature), color = "black") + 
        scale_y_discrete(labels = ClassA$drug) +
        scale_fill_manual(values = ARCHE_pal) +
        xlim(c(-1, 1)) + 
        theme_classic() + geom_vline(xintercept = 0) + 
        labs(y = "Drug", title = "", x = "Pearson's Correlation Coefficient", fill = "ARCHE")) 
    dev.off()
}

#' Plot Class B biomarker associations
#' 
plot_ClassB_biomarkersAssociations <- function(ClassB) {
    png("DrugResponse/results/figures/ClassB/biomarkersAssociations.png", width = 6, height = 5, res = 600, units = "in")
    print(ggplot(ClassB, aes(x = TE, y = rank)) +
        geom_col(aes(fill = signature), color = "black") + 
        scale_y_discrete(labels = ClassB$drug) +
        scale_fill_manual(values = ARCHE_pal) +
        xlim(c(-1, 1)) + 
        theme_classic() + geom_vline(xintercept = 0) + 
        labs(y = "Drug", title = "", x = "Transformed Effect Size", fill = "ARCHE")) 
    dev.off()

}

#' Forest-plot-like visualization for Class B
#'
#' Plot individual PC and meta-estimate for Class B biomarkers
#' 
plot_ClassB_forest <- function(toPlot) {

    png("DrugResponse/results/figures/ClassB/forest.png", width = 14, height = 6, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = estimate, y = pset)) + 
        geom_linerange(aes(xmin = lower, xmax = upper)) + 
        geom_vline(xintercept = 0) + 
        geom_vline(xintercept = c(-0.4, 0.4), linetype = "dotted") + 
        geom_point(data = toPlot, aes(shape = meta), size=6, fill = "#DB504A") +
        scale_shape_manual(values=c(23, 15), labels = c("Meta Estimate", "Pearson's\nCorrelation\nCoefficient")) + 
        scale_y_discrete(limits=rev) +
        guides(shape = guide_legend(title = NULL)) +
        facet_wrap(pairs ~ ., ncol = 7, nrow = 2) + 
        theme_classic() +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = "Pearson's Correlation Coefficient | Meta Estimate", y = "PSet"))
    dev.off()
}

#' Plot Class B associations as heatmap
#'
#' Annotated by ARCHE and BCa relevant drug
#' 
plot_ClassB_heatmap <- function(toPlot) {

    # format plot
    toPlot <- toPlot[order(toPlot$pairs),]
    toPlot$pset <- factor(toPlot$pset, levels = names(PSet_pal))
    toPlot$pairs <- factor(toPlot$pairs, levels = unique(toPlot$pairs))

    # arche bounds 
    n_pairs <- toPlot %>%
        group_by(signature) %>%
        summarise(n_pairs = n_distinct(pairs)) %>%
        mutate(cumulative = cumsum(n_pairs))
    bounds <- n_pairs$cumulative+0.5

    # main heatmap (pset~pair) PCC
    p1 <- ggplot(toPlot[toPlot$meta == FALSE,], aes(x = pairs, y = pset, fill = estimate)) + 
        geom_tile(color = 'black') +
        geom_text(data = subset(toPlot[toPlot$meta == FALSE,], FDR < 0.05),
                aes(label = "*"), 
                vjust = 0.75, size = 4) +
        geom_vline(xintercept = bounds[1:4], color = "gray") +
        scale_x_discrete(labels = function(pairs) gsub(".*_", "", pairs)) +
        scale_fill_gradient2("Pearson's\nCorrelation\nCoefficient", 
                            low = "#BC4749", 
                            high = "#689CB0",
                            mid = "#C2BBC9",
                            limits = c(-1, 1)) +
        theme_void() +
        theme(
            axis.text.y = element_text(size=9, hjust=1, vjust=0.5, margin = margin(r = 6)), 
            axis.text.x = element_text(size=9, angle=90, hjust=1, vjust=0.5, margin = margin(t = 3)),
            legend.title = element_text(size=9),
            axis.ticks = element_line(color = "gray", linewidth = 0.3),
            axis.ticks.length = unit(2, "pt")
            )
    
    # Meta-Estimate
    p2 <- ggplot(toPlot[toPlot$meta == TRUE,], aes(x = pairs, y = 1, fill = estimate)) + 
        geom_tile(color = "black") +
        theme_void() + 
        geom_vline(xintercept = bounds[1:4], color = "gray") +
        scale_fill_gradient2("Meta TE", 
                            low = "#BC4749", 
                            high = "#689CB0",
                            mid = "#C2BBC9",
                            limits = c(-1, 1)) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 5))) + 
        labs(y = "Meta TE")

    # BCa drug annotation
    p3 <- ggplot(toPlot, aes(x = pairs, y = 1, fill = ifelse(drug %in% bca_drugs, "A", "B"))) + 
        geom_tile(color = "gray") +
        theme_void() + 
        geom_vline(xintercept = bounds[1:4], color = "gray") +
        scale_fill_manual(values = c("#3E517A", "white")) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 2))) + 
        labs(y = "BCaDrug")
    
    # ARCHE annotation
    p4 <- ggplot(toPlot, aes(x = pairs, y = 1, fill = signature)) + 
        geom_tile(color = NA) +
        theme_void() + 
        geom_vline(xintercept = bounds[1:4], color = "gray") +
        scale_fill_manual("ARCHE", values = ARCHE_pal) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 5))) + 
        labs(y = "ARCHE ")

    # extract legends (just need the TE one, get others from Class A lol)
    l2 <- as_ggplot(get_legend(p2))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")
    p4 <- p4+theme(legend.position = "none")

    filename <- paste0("DrugResponse/results/figures/ClassB/heatmap.png")
    png(filename, width = 7, height = 3.5, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, p3, p4, l2, ncol = 10, nrow = 18,
        layout_matrix = rbind(c(4,4,4,4,4,4,4,4,4,4,4,4,4,5,5),
                              c(3,3,3,3,3,3,3,3,3,3,3,3,3,5,5),
                              c(2,2,2,2,2,2,2,2,2,2,2,2,2,5,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5),
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,NA,NA), 
                              c(1,1,1,1,1,1,1,1,1,1,1,1,1,NA,NA)))
    )
    dev.off()
}

#' Correlation plots of ARCHE and IC50
#' 
#' Helper function
#' 
plot_tdxd_corr <- function(combined, tdxd_PC, arche, ic50, save = FALSE) {

    # get corr
    pair <- paste(arche, ic50, sep = "_")
    corr <- round(tdxd_PC$pc[tdxd_PC$pairs == pair], 2)
    fdr <- round(tdxd_PC$FDR[tdxd_PC$pairs == pair], 2)

    p <- ggplot(combined, aes(x = .data[[arche]], y = .data[[ic50]], fill = subtype, color = subtype)) +
        geom_point(size = 3, shape = 21) +
        geom_smooth(method = "lm", se=F) + 
        scale_fill_manual(values = subtype_pal) +
        scale_color_manual(values = subtype_pal) +
        theme_classic() + 
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5, size = 10), 
            legend.key.size = unit(0.7, 'cm')
        ) +
        #xlim(-x, x) + ylim(-y, y) + 
        labs(title = paste0(arche, "\nPCC: ", corr, " (FDR: ", fdr, ")")) 

    if (save == TRUE) {
        filename <- paste0("DrugResponse/results/figures/TDXd/", arche, "_", ic50, ".png")
        png(filename, width = 5, height = 4, res = 600, units = "in")
        print(p)
        dev.off()
    } else {
        return(p)
    }
}

#' Correlation plots of ARCHE and IC50
#' 
#' Merge plots together
#' @param combined dataframe. Output of combined_data()
#' @param tdxd_PC dataframe. Output of computePC()
#' @param ic50 string. "Avg.IC50" or "Avg.IC50.Treps"
#' @param label string. For plot filename
#' 
plot_tdxd_all <- function(combined, tdxd_PC, ic50, label) {
    p1 <- plot_tdxd_corr(combined, tdxd_PC, "ARCHE1", ic50)
    p2 <- plot_tdxd_corr(combined, tdxd_PC, "ARCHE2", ic50)
    p3 <- plot_tdxd_corr(combined, tdxd_PC, "ARCHE3", ic50)
    p4 <- plot_tdxd_corr(combined, tdxd_PC, "ARCHE4", ic50)
    p5 <- plot_tdxd_corr(combined, tdxd_PC, "ARCHE5", ic50)
    p6 <- plot_tdxd_corr(combined, tdxd_PC, "ARCHE6", ic50)

    filename <- paste0("DrugResponse/results/figures/TDXd/", label, "_", ic50, ".png")
    png(filename, width = 10, height = 7, res = 600, units = "in")
    print(ggarrange(p1, p2, p3, p4, p5, p6, 
              nrow = 2, ncol = 3, 
              common.legend = TRUE, legend = "bottom"))
    dev.off()
}
