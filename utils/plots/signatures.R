#' ARCHE heatmap
#' 
#' Plot heatmap of ARCHE scores in TCGA BCa tumours
#' 
plot_ARCHE_heatmap <- function(mat) {

    omics_pal <- c("Yes" = random_blue, "No" = "white")

    # os
    p10 <- ggplot(mat, aes(x = rank, y = 1, fill = t_purity)) + 
        scale_fill_gradient(
            #name = "OS (Days)   ",
            #limits = c(130, 4300),
            name = "Tumour Purity",
            limits = c(0.4, 1),
            low = "#F8F4E3",
            high = "#7C3626",
            na.value = na_value
        ) +
        geom_tile(color = NA) +
        theme_void() + 
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    # stage
    p9 <- ggplot(mat, aes(x = rank, y = 1, fill = stage)) + 
        geom_tile(color = NA) +
        theme_void() + 
        scale_fill_manual(
            "Stage          ", 
            values = c(stage_pal), na.value = na_stage,
            labels = c("StageI", "StageII", "StageIII", "StageIV")) +
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    # age
    p8 <- ggplot(mat, aes(x = rank, y = 1, fill = age)) + 
        scale_fill_gradient(
            name = "Age           ",
            limits = c(30, 85),
            low = "#C1D4E6",
            high = "#13315C",
            na.value = na_value
        ) +
        geom_tile(color = NA) +
        theme_void() + 
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    # methylation
    p4 <- ggplot(mat, aes(x = rank, y = 1, fill = met)) + 
        geom_tile(color = "gray") +
        theme_void() + 
        scale_fill_manual("", values = omics_pal) +
        theme(axis.title.y = element_text(size=12), legend.position = "none") + 
        labs(y = "              ")

    # rna
    p5 <- ggplot(mat, aes(x = rank, y = 1, fill = rna)) + 
        geom_tile(color = "gray") +
        theme_void() + 
        scale_fill_manual("", values = omics_pal) +
        theme(axis.title.y = element_text(size=12), legend.position = "none") + 
        labs(y = "              ")

    # snv
    p6 <- ggplot(mat, aes(x = rank, y = 1, fill = snv)) + 
        geom_tile(color = "gray") +
        theme_void() + 
        scale_fill_manual("", values = omics_pal) +
        theme(axis.title.y = element_text(size=12), legend.position = "none") + 
        labs(y = "              ")
    
    # atac dup
    p7 <- ggplot(mat, aes(x = rank, y = 1, fill = dup)) + 
        geom_tile(color = "gray") +
        theme_void() + 
        scale_fill_manual("", values = omics_pal) +
        theme(axis.title.y = element_text(size=12), legend.position = "none") + 
        labs(y = "              ")

    # subtype
    p1 <- ggplot(mat, aes(x = rank, y = 1, fill = subtype)) + 
        geom_tile(color = NA) +
        theme_void() + 
        scale_fill_manual("BCa\nSubtype", values = subtype_pal,
                        labels = c("Basal", "Her2", "LuminalA", "LuminalB", "Normal", "NA")) +
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    # ARCHE
    p2 <- ggplot(mat, aes(x = rank, y = 1, fill = signature_assign)) + 
        geom_tile(color = NA) +
        theme_void() + 
        scale_fill_manual("Assigned\nARCHE", values = ARCHE_pal) +
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    # NMF output
    p3 <- ggplot(mat, aes(x = rank, y = ARCHE, fill = value)) + 
        geom_tile(color = NA) +
        scale_fill_gradientn("ARCHE        \nExpression\nScore", colours = brewer.pal(9, "Blues")) + 
        theme_void() +
        theme(
            axis.text.y = element_text(size=11, margin = margin(r = 2), hjust=1), 
            axis.title.x = element_text(size=12)
        ) + 
        labs(x = "Tumour Sample\n")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    l4 <- as_ggplot(get_legend(p8))
    l5 <- as_ggplot(get_legend(p9))
    l6 <- as_ggplot(get_legend(p10))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")
    p8 <- p8+theme(legend.position = "none")
    p9 <- p9+theme(legend.position = "none")
    p10 <- p10+theme(legend.position = "none")

    png("data/results/figures/1-Signatures/ARCHE_heatmap.png", width = 11, height = 6, res = 600, units = "in")
    print(
        grid.arrange(
            p1, p2, p3, l1, l2, l3, 
            p4, p5, p6, p7, p8, p9, 
            l4, l5, p10, l6, 
            ncol = 9, nrow = 13,
        layout_matrix = rbind(
                            c(10,10,10,10,10,10,10,NA,NA), 
                            c(9,9,9,9,9,9,9,NA,NA),
                            c(8,8,8,8,8,8,8,16,13), 
                            c(7,7,7,7,7,7,7,16,13),
                            c(11,11,11,11,11,11,11,16,13),
                            c(15,15,15,15,15,15,15,16,13),
                            c(12,12,12,12,12,12,12,16,13), 
                            c(1,1,1,1,1,1,1,4,14), 
                            c(2,2,2,2,2,2,2,4,14),
                            c(3,3,3,3,3,3,3,4,14), c(3,3,3,3,3,3,3,4,14),
                            c(3,3,3,3,3,3,3,4,14), c(3,3,3,3,3,3,3,4,14),
                            c(3,3,3,3,3,3,3,6,5), c(3,3,3,3,3,3,3,6,5),
                            c(3,3,3,3,3,3,3,6,5), c(3,3,3,3,3,3,3,6,5),
                            c(3,3,3,3,3,3,3,6,5), c(3,3,3,3,3,3,3,6,5)
        ))
    )
    dev.off()
}

#' Plot ARCHE window number
#' 
plot_ARCHE_peakInfo <- function(df, analysis) {

    df$ARCHE <- factor(df$ARCHE, levels = paste0("ARCHE", 6:1))
    filename <- paste0("data/results/figures/1-Signatures/ARCHE_", analysis, "_peakInfo.png")
    png(filename, width = 5, height = 5, res = 600, units = "in")
    print(
        ggplot(df, aes(x = num_windows, y = ARCHE)) +
            geom_bar(stat = "identity", fill = random_lightblue) +
            geom_text(data = df, aes(label = paste0(num_windows, " (", sum_peaks, "bp)"), x = 1.25e05)) +
            theme_classic() +
            labs(x = "Number of 500bp Windows")
    )
    dev.off()
}

#' Plot TCGA tumour staging by ARCHE
#' 
#' @param tnm string. Column label
#' 
plot_stage <- function(toPlot, tnm = "stage") {
    
    label <- switch(
        tnm,
        stage = "Stage",
        stageT = "Tumour Size",
        stageN = "Nodal Status",
        stageM = "Metastasis"
    )

    if (tnm == "stage") {
        p <- ggplot(toPlot, aes(x = signature_assign, fill = stage)) +
        geom_bar() +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(label, 
                    values = c(stage_pal), na.value = "#695F58",
                    labels = c("StageI", "StageII", "StageIII", "StageIV")) +
        theme_classic() + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(y = "Number of Tumours", x = "\nAssigned ARCHE")
    } else {
        p <- ggplot(toPlot, aes(x = signature_assign, fill = .data[[tnm]])) +
        geom_bar() +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(label, values = c(unname(stage_pal)), na.value = "#695F58") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(y = "Number of Tumours", x = "\nAssigned ARCHE")
    }

    png(paste0("data/results/figures/1-Signatures/ARCHE_", tnm, ".png"), width = 4, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot distribution of other clinical variabels by ARCHE
#' 
plot_clinical <- function(toPlot, var) {

    label <- switch(
        var,
        age = "Age",
        OS = "Overall Survival (Days)",
        t_purity = "Tumour Purity"
    )

    p <- ggplot(toPlot, aes(x = signature_assign, y = .data[[var]], fill = signature_assign)) +
    geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.5) +
    theme_classic() + 
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none") +
    scale_fill_manual(NULL, values = ARCHE_pal) +
    labs(y = label, x = "Assigned ARCHE")

    if (var == "t_purity") {
        p <- p + geom_signif(
            y_position = c(1.06, 1.02, 1, 0.91, 0.98),
            xmin = c(1, 2, 5, 3, 4), xmax = c(5, 5, 6, 5, 5),
            annotation = sig_codes,
            tip_length = 0.01,
            textsize = 4)
    }
    png(paste0("data/results/figures/1-Signatures/ARCHE_", var, ".png"), width = 4, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' ARCHE Upset plot
#' 
#' Plot Upset plot of overlapping peaks in ARCHEs
#' 
plot_ATAC_Upset <- function(m, label) {

    # set filter
    filter <- ifelse(label == "all", 2000000, 100000)
    m = m[comb_size(m) > filter]

    if (label == "20k") {
        w = 7
        h = 4
    } else {
        w = 11
        h = 5
    }

    filename <- paste0("data/results/figures/1-Signatures/ATAC_Upset_", label, ".png")
    png(filename, width = w, height = h, res = 600, units = "in")
    print(
        UpSet(m, set_order = c(paste0("ARCHE", 1:6)), comb_order = order(-comb_size(m)),
        bg_col = "gray", 
        right_annotation = upset_right_annotation(m))
    )
    dev.off()
}

#' Plot annotatePeak results
#' 
plot_annotatePeak <- function(toPlot, label) {

    filename <- paste0("data/results/figures/1-Signatures/annotatePeak_", label, ".png")
    png(filename, width = 7, height = 5, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(fill=Feature, y=Frequency, x=ARCHE)) + 
            geom_bar(position="fill", stat="identity", color = "black", size = 0.1) +
            scale_fill_manual(values = genfeat_pal) +
            theme_minimal() + 
            labs(y = "Percentage (%)")
        )
    dev.off()
}

#' Plot GREAT results
#' 
#' Plot lollipop graphs of GREAT enrcichment results
#' @param great dataframe. Output of runGREAT()
#' @param n_genes int. Number of annotated genes to filter by
#' @param label string. ARCHE label for filename
#' 
plot_GREAT <- function(great, n_genes, label) {

    # keep only top 30 with > n_genes
    great <- great[great$Total_Genes_Annotated > n_genes,]
    great <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
    toPlot <- great[1:30,]
    toPlot$name <- factor(toPlot$name, levels=rev(toPlot$name))

    p <- ggplot(toPlot, aes(x = name, y = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH)) +
        geom_point(aes(size = Total_Genes_Annotated), shape = 19) + 
        geom_segment(aes(x = name, xend = name, y = 0, yend = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH), size = 1) +
        coord_cartesian(clip = "off") + 
        coord_flip() +
        scale_alpha(range = c(1, 0.2)) +
        guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=19, size = 4))) +
        scale_colour_gradient(low = random_lightblue, high = "black") +
        theme_classic() + 
        theme(legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
              plot.title = element_text(hjust = 0.5, size = 12), 
              panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + 
        labs(x = "GO Term", y = "Fold Enrichment", size = "Genes", color = "Adjusted\nP-Value", title = label)

    filename <- paste0("data/results/figures/1-Signatures/GREAT/", label, "_ngene_", n_genes, ".png")
    png(filename, width = 8, height = 7, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot GREAT results with annotations
#' 
plot_GREAT_anno <- function(great, n_genes, label, anno) {

    # keep only top 20 with > n_genes
    great <- great[great$Total_Genes_Annotated > n_genes,]
    great <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
    toPlot <- great[1:20,]
    toPlot$name <- factor(toPlot$name, levels=rev(toPlot$name))

    # annotation bar
    anno$Label <- gsub("Developmental process", "Developmental\nprocess", anno$Label)
    anno <- anno[anno$ARCHE == label,]
    toPlot$anno <- anno$Label #[match(toPlot$ID, anno$GO)]

    p1 <- ggplot(toPlot, aes(x = name, y = 1, fill = anno)) +
        geom_tile(color = "black") +
        coord_flip() +
        scale_fill_manual(values = c(random_blue, "white")) +
        theme_void() +
        theme(
            axis.text.y = element_text(hjust = 1),
            axis.title.y = element_text(angle = 90),
            plot.margin=unit(c(1.03,-0.1,1.3,1), "cm")
        ) +
        labs(x = "GO Term", fill = "Related Pathways")

    # lollipop
    x <- max(toPlot$Hyper_Fold_Enrichment) + 0.25

    p2 <- ggplot(toPlot, aes(x = name, y = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH)) +
        geom_point(aes(size = Total_Genes_Annotated), shape = 19) + 
        geom_segment(aes(x = name, xend = name, y = 0, yend = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH), size = 1) +
        coord_cartesian(clip = "off") + 
        coord_flip() +
        scale_y_continuous(limits = c(0, x), expand = c(0,0)) +
        guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=19, size = 4))) +
        scale_colour_gradient(low = random_lightblue, high = "black") +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12), 
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.margin=unit(c(0.5,1,0.5,0), "cm")
        ) + 
        labs(y = "Fold Enrichment", size = "Number of Genes", color = "Adjusted P-Value", title = paste(label, "      "))

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")

    filename <- paste0("data/results/figures/1-Signatures/GREAT/", label, "_annotated_ngene_", n_genes, ".png")
    png(filename, width = 10, height = 7, res = 600, units = "in")
    print(grid.arrange(p1, p2, l1, l2, ncol = 7, nrow = 6,
        layout_matrix = rbind(c(1,1,1,1,2,2,NA),
                              c(1,1,1,1,2,2,3),
                              c(1,1,1,1,2,2,3),
                              c(1,1,1,1,2,2,4),
                              c(1,1,1,1,2,2,NA),
                              c(1,1,1,1,2,2,NA))))
    dev.off()

}

#' Plot HOMER motifs
#' 
plot_motifs <- function(arche, analysis) {
    
    label <- switch(
        analysis,
        known = paste(gsub("RCHE", "", arche), "K", sep = " - "),
        denovo = paste(gsub("RCHE", "", arche), "dn", sep = " - ")
    )

    file <- "data/results/data/1-Signatures/findMotifsGenome/findMotifsGenome.xlsx"
    motifs <- read_excel(file, sheet = label) |> as.data.frame()

    # get variables and format
    if (analysis == "known") {
        motifs$Name <- sub("\\/.*", "", motifs$Name)
    } else {
        motifs$Name <- sub("\\/.*", "", motifs$'Best Match/Details')
    }
    colnames(motifs) <- gsub(" Sequences with Motif", "", colnames(motifs))
    motifs <- motifs[,c("Name", "Rank", "P-value", "% of Targets", "% of Background")]
    motifs$Rank <- factor(motifs$Rank, levels = motifs$Rank)

    # save fold change results in another dataframe
    fc <- motifs
    fc$FoldChange <- motifs$'% of Targets' / motifs$'% of Background'

    motifs <- reshape2::melt(motifs, id.var = c("Name", "Rank", "P-value"))

    fc$logP <- log(fc$'P-value' + 1)

    # plot fc
    p <- ggplot(fc, aes(x = Rank, y = FoldChange, fill = .data[['% of Targets']])) + 
        #geom_point(aes(size = logP)) +
        geom_col(color = "black", linewidth = 0.2) +
        coord_flip() +
        scale_x_discrete(labels = setNames(motifs$Name, motifs$Rank)) +
        scale_fill_gradient(
            limits = c(0, 1),
            low = "#F8F4E3",
            high = "#7C3626",
            na.value = na_value
        ) +
        theme_minimal() +
        theme(
            legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 12), 
            legend.position = "none",
            axis.text.x = element_text(size = 10)
        ) +
        labs(
            x = "Transcription Factor (Family)",
            fill = "% of Target\nSequences\nwith Motif",
            title = paste(arche, analysis, "motifs"))

    filename <- paste0("data/results/figures/1-Signatures/motifs/", arche, "_", analysis, ".png")
    png(filename, width = 3.5, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot correlation plots of ARCHE scores and clinical variables
#' 
plot_corr_arche_cv <- function(mat, pheno) {

    arches <- paste0("ARCHE", 1:6)
    vars <- c("years_to_birth", "Tumor_purity", "overall_survival")

    for (arche in arches) {
        scores <- mat[mat$ARCHE == arche,]
        for (var in vars) {

            scores$p <- as.numeric(pheno[match(scores$variable, rownames(pheno)),][[var]])
            corr <- cor.test(scores$value, scores$p)
            pval <- round(corr$p.value, 4)
            pcc <- round(corr$estimate, 4)

            filename <- paste0("data/results/figures/1-Signatures/clinicalvars/", arche, "_", var, ".png")
            png(filename, width=5, height=4, units='in', res = 600, pointsize=80)
            print(ggplot(scores, aes(x = value, y = p)) +
                geom_point() +
                geom_smooth(method = "lm", se = FALSE) +
                theme_classic() + 
                theme(
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                    plot.title = element_text(size = 12)
                ) +
                labs(x = paste(arche, "Score"), y = var, title = paste0("PCC: ", pcc, "; pval: ", pval)))
            dev.off()
        }  
    }
}