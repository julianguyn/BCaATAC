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
get_all_drug_sen <- function() {

    # set up dataframe to store results
    all_drug_sen <- data.frame(matrix(nrow=0, ncol=3))
    colnames(all_drug_sen) <- c("Sample", "AAC", "PSet")
}

#' Get individual associations across PSets
#' 
#' @param pair string. In format "<Feature>_<Drug>"
#' @param signature_scores data.frame with Feature in rownames.
#' @param label string.
#' 
get_indiv_pcc <- function(pair, signature_scores, label, meta, pcc, plot = TRUE) {

    dir <- "4-DrugResponse/benchmarks/"

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
        x_label <- "Score"
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
        x_label <- "Score"
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
        x_label <- "Expression"
    }

    # print PCC
    for (pset in unique(all_drug_sen$PSet)) {
        subset <- all_drug_sen[all_drug_sen$PSet == pset,]
        if (nrow(subset) < 3) next
        pc <- cor.test(subset$Score, subset$AAC, method = "pearson", alternative = "two.sided")

        # update the global pcc dataframe
        pcc <<- rbind(pcc, data.frame(
            ARCHE_Drug = pair,
            Label = label,
            PSet = pset,
            PCC = round(pc$estimate, 4),
            pvalue = round(pc$p.value, 4)
        ))
    }

    # get x axis limits for plotting
    #x <- max(max(all_drug_sen$Score), abs(min(all_drug_sen$Score)))

    all_drug_sen$PSet <- factor(all_drug_sen$PSet, levels = names(PSet_pal))

    if (plot == TRUE) {
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
    p1 <- get_indiv_pcc(paste0(arche, "_", drug), arche_scores, "ARCHE", c_meta, pcc) + theme(legend.position = "none")
    p2 <- get_indiv_pcc(paste0(pam50, "_", drug), pam50_scores, "PAM50", c_meta, pcc) + theme(axis.title.y = element_blank())
    for (gene in unname(gene_list)) {
        get_indiv_pcc(paste0(gene, "_", drug), gene_df, "RNA", c_meta, pcc, plot = FALSE)
    }

    # plot ARCHE and PAM50
    message(length(unique(pcc$PSet)))
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
    names(cols) <- c(arche, drug, unname(gene_list))

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
    return(pcc)
}
