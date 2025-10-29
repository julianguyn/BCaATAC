#' Plot pheatmaps of mutational signatures
#' 
plot_pheatmap <- function(mat, label) {
    
    png(paste0("data/results/figures/2-MolecularSigAnalysis/molecularsig/pheatmap_", label, ".png"), width = 10, height = 10, res = 600, units = "in")
    print(
        pheatmap::pheatmap(mat = mat, 
                   cluster_rows = FALSE, 
                   main = paste("Cosine Similarity against", label, "Mutational Signatures"))
    )
    dev.off()
}

#' Plot heatmaps of mutational signatures
#' 

# helper function add ATAC signature to mutational sig dataframe
format_sig <- function(cosm_df) {
  cosm_df <- as.data.frame(cosm_df)
  cosm_df$Sample.Name <- rownames(cosm_df)
  cosm_df$ARCHE <- meta$ARCHE[match(cosm_df$Sample.Name, meta$Sample.Name)]
  cosm_df$rank <- meta$rank[match(cosm_df$Sample.Name, meta$Sample.Name)]
  cosm_df$Subtype <- meta$Subtype[match(cosm_df$Sample.Name, meta$Sample.Name)]
  cosm_df <- reshape2::melt(cosm_df, id = c("Sample.Name", "ARCHE", "Subtype", "rank"))
  cosm_df$rank <- factor(cosm_df$rank, levels = max(cosm_df$rank):1)
  return(cosm_df)
}

# plotting function
plot_heatmap_mutsig <- function(toPlot, type) {

    # add atac
    toPlot <- format_sig(toPlot)
  
    x_lab <- ifelse(type == "og30", "                     ", "                ")
    
    p1 <- ggplot(toPlot) + geom_tile(aes(x = variable, y = rank, fill = value), color = "gray") +
        scale_fill_distiller("Cosine Similarity", palette = "RdYlBu") +
        geom_hline(yintercept = c(16.5, 24.5, 40.5, 49.5, 58.5), linetype = "dotted", color = "black") +
        theme_void() + 
        theme(axis.text.x = element_text(size =8, angle = 90, vjust = 0.5, hjust=1), 
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(x = "\nMutational Signatures", y = "Tumour Sample\n")
    
    # signature annotation bar
    p2 <- ggplot(toPlot) + geom_tile(aes(x = 1, y = rank, fill = ARCHE)) +
        theme_void() +
        scale_fill_manual(values = ARCHE_pal) +
        theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(fill = "ARCHE          ", x = x_lab)
    
    # subtype annotation bar
    p3 <- ggplot(toPlot) + geom_tile(aes(x = 1, y = rank, fill = Subtype)) +
        theme_void() +
        scale_fill_manual(values = subtype_pal) +
        theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(fill = "Subtype      ", x = x_lab)
    
    
    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")
    
    png(paste0("data/results/figures/2-MolecularSigAnalysis/molecularsig/heatmap_", type, ".png"), width = 8, height = 6, res = 600, units = "in")
    print(grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 25, nrow = 6,
            layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4), 
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4), 
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5), 
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5),
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6),
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6))))
    dev.off()
}

#' Plot molecular sig - ARCHE correlations
#' 
plot_molecularsig_corr <- function(corr, label, type, data = "tumour") {
  
    corr$ATAC.Sig <- factor(corr$ATAC.Sig, levels = c(paste0("ARCHE", 1:6)))
    
    p <- ggplot(corr) + geom_tile(aes(x = Mol.Sig, y = ATAC.Sig, fill = Corr), color = "gray") +
        scale_fill_gradient2("Spearman\nCorrelation",
                            low = "#AF4C5B", high = "#046C9A", mid = "white",
                            midpoint = 0, limits = c(-1, 1)) +
        theme_void() + 
        theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1), 
            axis.text.y = element_text(size = 8, hjust = 0.95),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(x = paste0(label, "\n"), y = "\nARCHE") + coord_flip()

    # set width and height
    if (type == "og30") {
        w <- 4
        h <- 7
    } else if (type == "v3") {
        w <- 3
        h <- 8.5
    } else if (type == "hm") {
        w <- 5
        h <- 8
    } else if (type == "myc") {
        w <- 6
        h <- 8
    }
  
    png(paste0("data/results/figures/2-MolecularSigAnalysis/molecularsig/", type, "_", data, "_corr.png"), width = w, height = h, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot boxplots of molsig correlations with ARCHEs
#' 
plot_corr_boxplots <- function(corr_v3, corr_hm) {

    # format toPlot
    corr_v3$Signature <- "Mutational Signatures"
    corr_hm$Signature <- "Hallmark Gene Set"
    toPlot <- rbind(corr_v3[,c(1,3,4)], corr_hm[,c(1,3,4)])
    toPlot$Signature <- factor(toPlot$Signature, levels = c("Mutational Signatures", "Hallmark Gene Set"))

    # get min and max values 
    mm <- toPlot %>%
    group_by(ATAC.Sig, Signature) %>%
    filter(Corr == max(Corr) | Corr == min(Corr)) %>%
    ungroup()

    png("data/results/figures/2-MolecularSigAnalysis/corr_boxplots.png", width = 6, height = 4.5, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = ATAC.Sig, y = Corr, fill = Signature)) +
        geom_boxplot(color = "black") +
        theme_classic() + ylim(c(-1, 1)) + 
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
        scale_fill_manual(values = c("#077293", "#9AA5BD")) +
        geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
        labs(y = "Spearman Correlation", x = "") +
        geom_point(data = mm, 
             aes(x = ATAC.Sig, y = Corr, group = Signature), 
             color = "#CE3333", size = 2,
             position = position_dodge(width = 0.75),
             show.legend = FALSE)
    )
    dev.off()
}

