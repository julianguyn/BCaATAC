#' List of functions:
#' format_combinations()
#' compute_ci()
#' compute_pc()
#' compute_meta()

#' FDR correction and format dataframes 
#' 
#' Helper function for compute_ci() and compute_pc()
#' @param combinations dataframe. Intermediary from computeCI() and computePC().
#' @return A dataframe with FDR correction and additional formating.
#' 
format_combinations <- function(combinations, label) {
    
    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    #combinations$FDR <- p.adjust(combinations$pvalue, method = "BH", n = length(combinations$pvalue))
    combinations$FDR <- ave(combinations$pval, combinations$drug, FUN = function(p) p.adjust(p, method = "BH"))
    combinations$FDRsig <- ifelse(combinations$FDR < 0.05, TRUE, FALSE)

    # format dataframe for plotting (ordered by CI/PC already)
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$signature, "_", combinations$drug)
    combinations$pset <- c(rep(label, nrow(combinations)))

    return(combinations)
}

#' Compute Pearson's correlation
#'
#' Computes PC between ARCHE scores and drug AAC.
#' @param signature_scores dataframe. ARCHE scores from get_scores().
#' @param sensitivity_data dataframe. Drug sensitivity from get_drugsen().
#' @param label string. PSet name
#' @return A dataframe of CI and metrics.
#' 
compute_pc <- function(signature_scores, sensitivity_data, label) {
    
    # get PSet-specific cells
    to_keep <- intersect(colnames(signature_scores),colnames(sensitivity_data))
    signature_scores <- signature_scores[,match(to_keep, colnames(signature_scores))]
    sensitivity_data <- sensitivity_data[,match(to_keep, colnames(sensitivity_data))]

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("signature", "drug", "pc", "pvalue", "n", "upper", "lower")
    combinations$signature <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){
        AAC <- as.numeric(sensitivity_data[combinations$drug[i],])
        if (length(AAC[!is.na(AAC)]) > 3) {
            pc <- cor.test(AAC, # drug AAC
                    as.numeric(unlist(signature_scores[combinations$signature[i],])), # ARCHE score
                    method = 'pearson', alternative = 'two.sided')

            combinations$pvalue[i] <- pc$p.value
            combinations$pc[i] <- pc$estimate
            combinations$n[i] <- length(AAC)
            combinations$upper[i] <- pc$conf.int[1]
            combinations$lower[i] <- pc$conf.int[2]
        }
    }

    # format dataframe and FDR correction
    combinations <- combinations[order(combinations$pc),]
    combinations <- format_combinations(combinations, label)
    
    return(combinations)
}

#' Perform meta analysis
#'
#' Compute meta estimates for signature drug response associations in >=3 PSets.
#' @param df dataframe. Dataframe from computePC()
#' @return A dataframe of meta analysis results.
#' 
compute_meta <- function(df) {
    
    # keep only signature-drug pairs that are present in at least 3 PSets
    df <- df[which(df$pairs %in% names(table(df$pairs)[table(df$pairs) > 2])),]

    # data frame to hold meta estimates
    estimates <- as.data.frame(matrix(data = NA, nrow = length(unique(df$pairs)), ncol = 10))
    colnames(estimates) <- c("signature","drug", "pair", "TE", "seTE", "upper", "lower", "pval", "pval.Q", "I2")

    # perform meta-analysis
    for (i in 1:length(unique(df$pair))) {
        
        pair <- unique(df$pairs)[i]
        tmp <- df[which(df$pair == pair),]
        meta <- metacor(pc, n, data = tmp, method.tau = "DL", studlab = tmp$pset)
        
        estimates$signature[i] <- tmp$signature[1]
        estimates$drug[i] <- tmp$drug[1]
        estimates$pair[i] <- pair
        estimates$TE[i] <- meta$TE.random
        estimates$seTE[i] <- meta$seTE.random
        estimates$upper[i] <- meta$upper.random
        estimates$lower[i] <- meta$lower.random
        estimates$pval[i] <- meta$pval.random
        estimates$pval.Q[i] <- meta$pval.Q
        estimates$I2[i] <- meta$I2

        if (abs(meta$TE.random) > 0.4) {
            # plot forest plot
            fileName = paste0("data/results/figures/4-DrugResponse/ClassB/meta/",pair,".png")
            
            png(fileName, width = 10, height = 4, res = 600, units = "in")
            title <- pair
            forest(meta,
                leftcols = c("studlab", "TE", "seTE", "lower", "upper", "pval"),
                leftlabs = c(title, "Effect", "SE", "95% CI \n Lower", "95% CI \n Upper", "P value"),
                xlab = "effect estimate", lab.e = "Intervention", sortvar = TE, smlab = " ", text.random = "Random effect", 
                print.I2.ci = FALSE, print.Q = TRUE, print.pval.Q = TRUE, digits.sd = 2, print.I2 = TRUE, print.tau2 = TRUE,
                text.random.w = TRUE, colgap.forest.left = "0.5cm", layout = "RevMan5", test.overall.random = TRUE,
                test.overall.common = TRUE, xlim = "symmetric", col.square = "grey70", col.inside = "grey70", col.square.lines = "grey30", 
                col.diamond.random = "#526863", col.diamond.common  = "#BD6B73", ff.xlab = "bold", fontsize = 11, fs.heading = 11.5,
                squaresize = 0.55, scientific.pval = TRUE, lty.random = NULL, lty.fixed  = NULL)
            dev.off()
        }
    }

    estimates$FDR <- ave(estimates$pval, estimates$drug, FUN = function(p) p.adjust(p, method = "BH"))
    
    return(estimates)
}

#' Helper function to compute ttest of drug response
#' 
#' @param sen sensitivity dataframe
#' @param pset string label
#' @param tt_res dataframe to store results
#' @param l1_cells string cells in group1
#' @param l2_cells string cells in group2
#' @param lab1 string ARCHE label for group1
#' @param lab2 string ARCHE label for group2
#' 
diff_drug_response <- function(sen, pset, tt_res, l1_cells, l2_cells, lab1, lab2) {
    sen1 <- sen[,colnames(sen) %in% l1_cells]
    sen2 <- sen[,colnames(sen) %in% l2_cells]

    for (drug in rownames(sen1)) {
        sen1_AAC <- na.omit(as.numeric(sen1[drug, ]))
        sen2_AAC <- na.omit(as.numeric(sen2[drug, ]))

        if (length(sen1_AAC) >= 2 && length(sen2_AAC) >= 2) { 
            tt <- t.test(sen1[drug,], sen2[drug,])
            res <- data.frame(
                Drug = drug,
                PSet = pset,
                pval = tt$p.value,
                statistic = tt$statistic,
                group1_mean = tt$estimate[1],
                group2_mean = tt$estimate[2],
                sensitive = ifelse(tt$estimate[1] > tt$estimate[2], lab1, lab2)
            )
            tt_res <- rbind(tt_res, res)
        }
    }
    return(tt_res)
}

#' Helper function to plot boxplots
plot_diff_boxplots <- function(toPlot, label) {
    toPlot$Drug <- sub(":", ":\n", sub(" \\(.*\\)", "", toPlot$Drug))
    toPlot$Drug <- paste0(sub("_", "\n(", toPlot$Drug), ")")
    p <- ggplot() +
        geom_boxplot(
            data = toPlot,
            aes(x = Drug, y = AAC, fill = ARCHE),
            outlier.shape = NA) +
        scale_fill_manual(values = ARCHE_pal) +
        new_scale_fill() +
        geom_jitter(
            data = toPlot, aes(x = Drug, y = AAC, group = ARCHE, fill = Subtype),
            position = position_jitterdodge(jitter.width = 0.5),
            shape = 21, size = 2) +
        scale_fill_manual(values = subtype_pal) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
        ) +
        ylim(c(0, 0.7)) +
        labs(title = label)
    return(p)
}

#' Compile drug response difference to compare between ARCHEs
#' 
compile_diff <- function(arche_scores, lab1, lab2, drug_pal, label) {

    tt <- as.data.frame(t(arche_scores[c(lab1, lab2),]))
    l1_cells <- rownames(tt[tt[[lab1]] > quantile(tt[[lab1]])["75%"],])
    l2_cells <- rownames(tt[tt[[lab2]] > quantile(tt[[lab2]])["75%"],])

    to_rm <- intersect(l1_cells, l2_cells) # two cells: CAL51 and SUM149PT, 10 unique in each
    cat("Removing:\n", to_rm, "\n")
    l1_cells <- l1_cells[-which(l1_cells %in% to_rm)]
    l2_cells <- l2_cells[-which(l2_cells %in% to_rm)]
    cat("Number of cells from", lab1, ":", length(l1_cells), "\n")
    cat("Number of cells from", lab2, ":", length(l2_cells), "\n")

    # compile all drug response
    tt_res <- data.frame(matrix(nrow = 0, ncol = 7))

    tt_res <- diff_drug_response(ubr1_sen, "UBR1", tt_res, l1_cells, l2_cells, lab1, lab2)
    tt_res <- diff_drug_response(ubr2_sen, "UBR2", tt_res, l1_cells, l2_cells, lab1, lab2)
    tt_res <- diff_drug_response(gray_sen, "GRAY", tt_res, l1_cells, l2_cells, lab1, lab2)
    tt_res <- diff_drug_response(gcsi_sen, "gCSI", tt_res, l1_cells, l2_cells, lab1, lab2)
    tt_res <- diff_drug_response(gdsc_sen, "GDSC2", tt_res, l1_cells, l2_cells, lab1, lab2)
    tt_res <- diff_drug_response(ccle_sen, "CCLE", tt_res, l1_cells, l2_cells, lab1, lab2)
    tt_res <- diff_drug_response(ctrp_sen, "CTRP", tt_res, l1_cells, l2_cells, lab1, lab2)


    sig <- tt_res[tt_res$pval < 0.1,]
    sig$Diff <- abs(sig$group1_mean - sig$group2_mean)
    sig <- sig[sig$Diff > 0.1,]

    toPlot <- data.frame(
        Drug = rep(paste0(sig$Drug, "_", sig$PSet), 2),
        ARCHE = c(rep(lab1, nrow(sig)), rep(lab2, nrow(sig))),
        Mean = c(sig$group1_mean, sig$group2_mean)
    )
    toPlot$Sensitive <- sig$sensitive[match(sub("_.*", "", toPlot$Drug), sig$Drug)]

    toPlot$Label <- ifelse(
        paste0(sub("_.*", "", toPlot$Drug) %in% names(drug_pal)),
        paste0(sub("_.*", "", toPlot$Drug)), "")
    toPlot$Label <- factor(toPlot$Label, levels = c(names(drug_pal), ""))

    # plot dot plots
    p <- ggplot(toPlot, aes(x = ARCHE, y = Mean, group = Drug)) +
        geom_line(alpha = 0.5, color = "gray") +
        geom_point(size = 4, color = "gray") +
        facet_wrap(~Sensitive, ncol = 2) +
        new_scale_fill() +
        geom_line(
            data = toPlot[toPlot$Label != "",],
            aes(x = ARCHE, y = Mean, group = Drug),
            alpha = 0.8, color = "black") +
        geom_point(
            data = toPlot[toPlot$Label != "",],
            aes(x = ARCHE, y = Mean, group = Drug, fill = Label),
            size = 4, shape = 21) +
        scale_fill_manual(
            values = drug_pal,
            labels = sub(" ", "\n", names(drug_pal))) +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = "white"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
            axis.title.y = element_text(margin = margin(r = 10))
        ) +
        labs(y = "Average AAC")
    filename <- paste0("data/results/figures/4-DrugResponse/comp/", label, "_dots.png")
    ggsave(filename, p, w=5, h=4)
    

    # get all drug response
    all_response <- data.frame(matrix(nrow=0, ncol=5))
    for (drug in names(drug_pal)) {
        response <- compile_AAC(drug)
        response <- response[response$Sample %in% c(l1_cells, l2_cells),]
        response$ARCHE <- ifelse(response$Sample %in% l1_cells, lab1, lab2)
        response$Drug <- drug
        all_response <- rbind(all_response, response)
    }
    all_response$Drug <- paste0(all_response$Drug, "_", all_response$PSet)
    all_response <- all_response[all_response$Drug %in% toPlot$Drug,]
    all_response$Sensitive <- toPlot$Sensitive[match(all_response$Drug, toPlot$Drug)]
    all_response$Subtype <- c_meta$subtype[match(all_response$Sample, c_meta$sampleid)]


    p1 <- plot_diff_boxplots(toPlot <- all_response[grep("Clofarabine|Cytarabine|Gemcitabine", all_response$Drug),], "Antimetabolites (Single)")
    p2 <- plot_diff_boxplots(toPlot <- all_response[grep("mol", all_response$Drug),], "Antimetabolites (Combo)")
    p3 <- plot_diff_boxplots(toPlot <- all_response[grep("SN-38|Topotecan", all_response$Drug),], "TOPI inhibitor") + theme(axis.title.y = element_blank())
    p4 <- plot_diff_boxplots(toPlot <- all_response[grep("Dasatinib", all_response$Drug),], "Dasatinib") + theme(axis.title.y = element_blank())

    p_upper <- p1 + p3 + plot_layout(widths = c(1.5, 1), guides = "collect") &
        theme(legend.position = "right")
    p_lower <- p2 + p4 + plot_layout(widths = c(1.5, 1), guides = "collect") &
        theme(legend.position = "right")

    p <- wrap_elements(p_upper) / wrap_elements(p_lower)
    filename <- paste0("data/results/figures/4-DrugResponse/comp/", label, "_boxplots.png")
    ggsave(filename, p, w=6.5, h=5)

}
