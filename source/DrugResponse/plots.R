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
            labels = c("Positive\n(Concordant)\nAssociation",
                    "Negative\n(Discordant)\nAssociation")) +
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

#' Plot Class A associations in >1 PSet as heatmap
#'
#' Annotated by ARCHE and BCa relevant drug
#' 
plot_ClassA_heatmapMulti <- function(toPlot) {

    # subset signature and format plot
    toPlot <- toPlot[order(toPlot$pairs),]
    toPlot$pset <- factor(toPlot$pset, levels = names(PSet_pal))
    toPlot$pairs <- factor(toPlot$pairs, levels = rev(unique(toPlot$pairs)))

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
                            low = binary_pal[2], 
                            high = binary_pal[1],
                            limits = c(-1, 1)) +
        theme_void() +
        theme(
            axis.text.y = element_text(size=9, hjust=1, vjust=0.5, margin = margin(r = 6)), 
            axis.text.x = element_text(size=9, angle=90, hjust=1, vjust=0.5, margin = margin(t = 3)),
            legend.title = element_text(size=9),
            axis.ticks = element_line(color = "gray", linewidth = 0.3),
            axis.ticks.length = unit(2, "pt")
            )

    # ARCHE annotation
    p2 <- ggplot(toPlot, aes(x = pairs, y = 1, fill = signature)) + 
        geom_tile(color = NA) +
        theme_void() + 
        geom_vline(xintercept = bounds[1:5], color = "gray") +
        scale_fill_manual("ARCHE", values = ARCHE_pal) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 5))) + 
        labs(y = "ARCHE ")

    # BCa drug annotation
    p3 <- ggplot(toPlot, aes(x = pairs, y = 1, fill = ifelse(drug %in% bca_drugs, "A", "B"))) + 
        geom_tile(color = NA) +
        theme_void() + 
        geom_vline(xintercept = bounds[1:5], color = "gray") +
        scale_fill_manual(values = c("#3E517A", "white")) +
        theme(axis.title.y = element_text(size = 9, hjust=1, margin = margin(r = 2))) + 
        labs(y = "BCaDrug")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")

    png("DrugResponse/results/figures/ClassA/heatmapMulti.png", width = 18, height = 3.5, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, p3, l2, l1, ncol = 10, nrow = 17,
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
    png("DrugResponse/results/figures/ClassA/biomarkersAssociations.png", width = 6, height = 5, res = 600, units = "in")
    print(ggplot(ClassA, aes(x = pc, y = rank)) +
        geom_col(aes(fill = signature), color = "black") + 
        scale_y_discrete(labels = ClassA$drug) +
        scale_fill_manual(values = ARCHE_pal) +
        xlim(c(0, 1)) + 
        theme_classic() + geom_vline(xintercept = 0) + 
        labs(y = "Drug", title = "", x = "Pearson's Correlation Coefficient", fill = "ARCHE")) 
    dev.off()

}

#' Plot Class A associations across PSets
#'
#' Colour by if == ClassA biomarker
#' 
plot_ClassA_associationsAcrossPSets <- function(toPlot) {

    png("DrugResponse/results/figures/ClassA/associationsAcrossPSets.png", width = 17, height = 5, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = pset, y = pc, fill = ifelse((FDR<0.05 & abs(pc) > 0.6), "Y", "N"))) +
        geom_bar(stat="identity", color = "black") +
        geom_text(data = subset(toPlot, FDR < 0.05),
            aes(label = "*", y = pc), 
            vjust = 0, size = 6) +
        facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
        scale_fill_manual(values = c("#8BA6A9", "#6F5E5C"), breaks = c("Y", "N"),
            labels = c("ClassA:\nabs(PCC>=0.65)\n& FDR<0.05",
                    "Not ClassA:\nabs(PCC<0.65)\nor FDR>0.05")) +
        geom_hline(yintercept = c(0.6), linetype = "dotted") +
        ylim(c(0,1)) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.spacing = unit(0, "lines")
        ) +
        labs(y = "Pearson's Correlation Coefficient", fill = "Biomarker\nAssociation"))
    dev.off()
}

#' Plot Class B biomarker associations
#' 
plot_ClassB_biomarkersAssociations <- function(ClassB) {

    toPlot <- ClassB[which(ClassB$FDR < 0.05 & abs(ClassB$pc) > 0.4),]

    png("DrugResponse/results/figures/ClassB/biomarkersAssociations.png", width = 17, height = 5, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = pset, y = pc, fill = signature)) +
        geom_bar(stat="identity", color = "black") +
        facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
        scale_fill_manual(values = ARCHE_pal) +
        geom_hline(yintercept = c(0.4, -0.4), linetype = "dotted") +
        ylim(c(0,1)) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.spacing = unit(0, "lines")
        ) +
        labs(y = "Pearson's Correlation Coefficient", fill = "Biomarker\nAssociation"))
    dev.off()
}


#' Plot Class B associations across PSets
#'
#' Colour by if == ClassB biomarker
#' 
plot_ClassB_associationsAcrossPSets <- function(toPlot) {

    png("DrugResponse/results/figures/ClassB/associationsAcrossPSets.png", width = 17, height = 5, res = 600, units = "in")
    print(ggplot(ClassB, aes(x = pset, y = pc, fill = ifelse((FDR<0.05 & abs(pc) > 0.4), "Y", "N"))) +
        geom_bar(stat="identity", color = "black") +
        geom_text(data = subset(ClassB, FDR < 0.05),
            aes(label = "*", y = pc), 
            vjust = 0, size = 6) +
        facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
        scale_fill_manual(values = c("#8BA6A9", "#6F5E5C"), breaks = c("Y", "N"),
            labels = c("ClassB:\nabs(PCC>=0.4)\n& FDR<0.05",
                    "Not ClassB:\nabs(PCC<0.4)\nor FDR>0.05")) +
        geom_hline(yintercept = c(0.4, -0.4), linetype = "dotted") +
        geom_hline(yintercept = 0) +
        ylim(c(-0.5,1)) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.spacing = unit(0, "lines")
        ) +
        labs(y = "Pearson's Correlation Coefficient", fill = "Biomarker\nAssociation"))
    dev.off()
}

#' Plot Class C biomarker associations
#' 
plot_ClassC_biomarkersAssociations <- function(ClassC) {
    png("DrugResponse/results/figures/ClassC/biomarkersAssociations.png", width = 6, height = 5, res = 600, units = "in")
    print(ggplot(ClassC, aes(x = TE, y = rank)) +
        geom_col(aes(fill = signature), color = "black") + 
        scale_y_discrete(labels = ClassA$drug) +
        scale_fill_manual(values = ARCHE_pal) +
        xlim(c(-1, 1)) + 
        theme_classic() + geom_vline(xintercept = 0) + 
        labs(y = "Drug", title = "", x = "Transformed Effect Size", fill = "ARCHE")) 
    dev.off()

}

#' Forest-plot-like visualization for Class C
#'
#' Plot individual PC and meta-estimate for Class C biomarkers
#' @param df dataframe. Dataframe from computePC()
#' @param ClassC dataframe. Subsetted dataframe from run_meta() of ClassC pairs
#' 
plot_ClassC_forest <- function(PC_res, ClassC) {

    # only plot classC
    keep <- PC_res[which(PC_res$pairs %in% ClassC$pair),]
    
    # combine PC_res and meta results
    toPlot <- data.frame(
        signature = c(keep$signature, ClassC$signature),
        drug = c(keep$drug, ClassC$drug),
        pairs = c(keep$pairs, ClassC$pair),
        estimate = c(keep$pc, ClassC$TE),
        upper = c(keep$upper, ClassC$upper),
        lower = c(keep$lower, ClassC$lower),
        FDR = c(keep$FDR, ClassC$FDR),
        pset = c(keep$pset, rep("Meta Estimate", nrow(ClassC)))
    )
    toPlot$meta <- factor(ifelse(toPlot$pset == "Meta Estimate", TRUE, FALSE), levels = c(TRUE, FALSE))
    toPlot$pset <- factor(toPlot$pset, levels = c(unique(PC_res$pset), "Meta Estimate"))

    # plot
    png("DrugResponse/results/figures/ClassC/forest.png", width = 14, height = 6, res = 600, units = "in")
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
