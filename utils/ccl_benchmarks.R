#' Helper function to compile pam50 scores
format_pam50 <- function(pam50, label) {
    pam50$Label <- label
    pam50$Sample <- rownames(pam50)
    rownames(pam50) <- NULL
    return(pam50)
}

#' Helper function to compile gene expression
format_rna <- function(rna, label, gene_list) {
    rna <- rna[names(gene_list),]
    rownames(rna) <- gene_list[rownames(rna)]
    rna <- as.data.frame(t(rna))
    rna$PSet <- label
    rna$Sample <- rownames(rna)
    rownames(rna) <- NULL
    return(rna)
}

#' Subfunction to compile all_drug_sen
get_all_drug_sen <- function(pair, signature_scores, label, meta) {

    message(paste("Getting drug_sen for", pair))
    feature <- gsub("_.*", "", pair)
    drug <- gsub(paste0(feature, "_"), "", pair)

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
    all_drug_sen$Score <- NA
    if (label == "ARCHE") {
        signature_scores <- signature_scores[rownames(signature_scores) == feature,]
        all_drug_sen <- all_drug_sen[all_drug_sen$Sample %in% colnames(signature_scores),]
        for (i in 1:nrow(all_drug_sen)) { 
            all_drug_sen$Score[i] <-  signature_scores[,colnames(signature_scores) == all_drug_sen$Sample[i]]
        }
        all_drug_sen$Subtype <- meta$subtype[match(all_drug_sen$Sample, meta$sampleid)]
    } else if (label == "PAM50") {
        all_drug_sen$Subtype <- NA
        all_drug_sen <- all_drug_sen[all_drug_sen$Sample %in% signature_scores$Sample,]
        for (i in 1:nrow(all_drug_sen)) {
            pset <- all_drug_sen$PSet[i]
            cell <- all_drug_sen$Sample[i]
            tt <- signature_scores[signature_scores$Sample == cell & signature_scores$Label == pset,]
            if (nrow(tt) > 0) {
                all_drug_sen$Score[i] <- tt[[feature]]
                all_drug_sen$Subtype[i] <- tt$Subtype
            }
        }
        all_drug_sen <- all_drug_sen[!is.na(all_drug_sen$Score), ]
    } else if (label == "RNA") {
        all_drug_sen <- all_drug_sen[all_drug_sen$Sample %in% signature_scores$Sample,]
        for (i in 1:nrow(all_drug_sen)) {
            pset <- all_drug_sen$PSet[i]
            cell <- all_drug_sen$Sample[i]
            tt <- signature_scores[signature_scores$Sample == cell & signature_scores$PSet == pset,]
            if (nrow(tt) > 0) {
                all_drug_sen$Score[i] <- tt[[feature]]
            }
        }
        all_drug_sen <- all_drug_sen[!is.na(all_drug_sen$Score), ]
        all_drug_sen$Subtype <- meta$subtype[match(all_drug_sen$Sample, meta$sampleid)]
    }
    return(all_drug_sen)
}

#' Subfunction to plot faceted scatter plots
plot_scatter <- function(pair, all_drug_sen, label) {

    message(paste("Plotting scatter for", pair))

    feature <- gsub("_.*", "", pair)
    drug <- gsub(paste0(feature, "_"), "", pair)
    x_label <- ifelse(label == "RNA", "Expression", "Score")

    all_drug_sen$PSet <- factor(all_drug_sen$PSet, levels = names(PSet_pal))
    p <- ggplot(all_drug_sen, aes(x = Score, y = AAC, fill = Subtype)) + 
        geom_point(size = 2.5, shape = 21) + 
        facet_wrap(~PSet, nrow = length(unique(all_drug_sen$PSet))) +
        geom_smooth(method = "lm", se=TRUE, color = "black", aes(group = 1), show.legend = FALSE) + 
        scale_fill_manual(values = subtype_pal) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        theme_bw() + 
        theme(
            strip.background = element_rect(fill = "white"),
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), 
            legend.key.size = unit(0.5, 'cm'),
            axis.title.y = element_text(size = 11, margin = margin(r = 10))
        ) +
        #xlim(-x, x) + 
        ylim(0, 1) + 
        labs(
            x = paste(feature, x_label), 
            y = paste(drug, "Response (AAC)"),
            title = label
        )
    return(p)
}

#' Get individual associations across PSets
#' 
#' @param pair string. In format "<Feature>_<Drug>"
#' @param all_drug_sen data.frame from get_all_drug_sen
#' @param label string.
#' 
get_pcc <- function(pair, all_drug_sen, pcc) {

    message(paste("Getting PCC for", pair))

    for (pset in unique(all_drug_sen$PSet)) {
        subset <- all_drug_sen[all_drug_sen$PSet == pset,]
        if (nrow(subset) < 3) next
        pc <- cor.test(subset$Score, subset$AAC, method = "pearson", alternative = "two.sided")

        pcc <- rbind(pcc, data.frame(
            Feature_Drug = pair,
            PSet = pset,
            PCC = round(pc$estimate, 4),
            pvalue = round(pc$p.value, 4)
        ))
    }
    return(pcc)

}

#' Plot all associations for given drug
#' 
#' @param arche string. ARCHE name
#' @param pam50 string. PAM50 subtype name
#' @param drug string. Drug
#' @param arche_scores dataframe. ARCHE scores
#' @param gene_list vector of genes
#' 
plot_associations <- function(arche, pam50, drug, arche_scores, gene_list) {

    # compile gene expression
    gene_df <- rbind(
        format_rna(ubr1, "UBR1", gene_list), format_rna(ubr2, "UBR2", gene_list),
        format_rna(gray, "GRAY", gene_list), format_rna(gcsi, "gCSI", gene_list),
        format_rna(gdsc, "GDSC2", gene_list), format_rna(ccle, "CCLE", gene_list),
        format_rna(ccle, "CTRP", gene_list)
    )

    # create dataframe to store results
    pcc <- data.frame(matrix(nrow=0, ncol=5))
    colnames(pcc) <- c("Feature_Drug", "Label", "PSet", "PCC", "pvalue")

    # get individual associations
    arche_sen <- get_all_drug_sen(paste0(arche, "_", drug), arche_scores, "ARCHE", c_meta)
    pam50_sen <- get_all_drug_sen(paste0(pam50, "_", drug), pam50_scores, "PAM50", c_meta)

    pcc <- get_pcc(paste0(arche, "_", drug), arche_sen, pcc)
    pcc <- get_pcc(paste0(pam50, "_", drug), pam50_sen, pcc)
    for (gene in unname(gene_list)) {
        pcc <- get_pcc(paste0(gene, "_", drug), get_all_drug_sen(paste0(gene, "_", drug), gene_df, "RNA", c_meta), pcc)
    }
    message(nrow(pcc))
    message(ncol(pcc))

    # plot ARCHE and PAM50
    p1 <- plot_scatter(paste0(arche, "_", drug), arche_sen, "ARCHE") + theme(legend.position = "none")
    p2 <- plot_scatter(paste0(pam50, "_", drug), pam50_sen, "PAM50") + theme(axis.title.y = element_blank())
    p <- p1 + p2
    filename <- paste0("data/results/figures/4-DrugResponse/benchmarks/", drug, "/arche_pam50.png")
    ggsave(filename, p, width = 5.5, height = length(unique(pcc$PSet))+3)

    # format toPlot
    pcc$Feature <- sub("_.*", "", pcc$Feature_Drug)
    pcc$Feature <- factor(pcc$Feature, levels = unique(pcc$Feature))
    pcc$sig <- ifelse(pcc$pvalue < 0.01, "P-Value < 0.05", "P-Value >= 0.05")
    pcc$PSet <- factor(pcc$PSet, levels = names(PSet_pal))

    # set colour palette
    cols <- c("#3F3244", "#208AAE", rep("#BFC3BA", length(gene_list)))
    names(cols) <- c(arche, pam50, unname(gene_list))

    # plot all markers
    p <- ggplot(pcc, aes(x = PSet, y = PCC, fill = Feature, size = -log(pvalue))) +
        geom_point(alpha = 0.7, shape = 21) +
        scale_size(range = c(1, 12)) +
        scale_fill_manual(values = cols) +
        guides(fill = guide_legend(override.aes = list(size = 3))) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            legend.key.size = unit(0.4, 'cm'),
            legend.title = element_text(size = 10)
        )
    filename <- paste0("data/results/figures/4-DrugResponse/benchmarks/", drug, "/all_compiled.png")
    ggsave(filename, p, width = 5.5, height = 3)
    rownames(pcc) <- NULL
    return(pcc)
}
