#' Plot bvals per ARCHE
#' 
plot_bval_arches <- function(toPlot) {

    # categorize methylation status
    toPlot <- toPlot %>%
        mutate(
            status = case_when(
                value <  0.2 ~ "Hypomethylated",
                value > 0.8 ~ "Hypermethylated",
                between(value, 0.2, 0.8) ~ "Intermediate",
                TRUE ~ NA_character_
            )
        )

    # get proportions by methylation status
    toPlot <- toPlot %>%
        filter(!is.na(status)) %>%
        group_by(ARCHE, label, status) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(ARCHE, label) %>%
        mutate(prop = n / sum(n))
    toPlot$status <- factor(toPlot$status, levels = c("Hypomethylated", "Intermediate", "Hypermethylated"))

    # plot
    filename <- "data/results/figures/2-MolecularSigAnalysis/methylation/bval_arches.png"
    png(filename, width = 6, height = 4, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = label, y = prop, fill = status)) +
            geom_bar(stat = "identity", position = "stack") +
            facet_wrap(~ ARCHE) +
            geom_text(aes(label = paste0(round(prop * 100), "%")), position = position_fill(vjust = 0.5)) +
            scale_fill_manual(values = c(random_blue, "gray", random_lightblue)) +
            theme_classic() +
            theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
            labs(x = "", y = "Proportion of CpG sites", fill = "Methylation Status")
    )
    dev.off()
}

#' Plot DMPs
#' 
plot_dmps_arches <- function(toPlot, counts) {
    filename <- "data/results/figures/2-MolecularSigAnalysis/methylation/dmps_arches.png"
    png(filename, width = 6, height = 4, res = 600, units = "in")
    print(
        ggplot() +
            geom_col(data = toPlot, aes(x = label, y = n, fill = anno)) +
            geom_text(data = counts, aes(x = label, y = n, label = n), vjust = -0.5) +
            facet_wrap(~ARCHE) +
            scale_y_continuous(limits = c(0, 117), expand = c(0, 0)) +
            scale_fill_manual("CpG Island\nMapping", values = dmp_pal) +
            theme_classic() + 
            theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
            labs(y = "Count", x = "")
    )
    dev.off()
}

#' Plot mean diff of all DMRs for an ARCHE
#' 
plot_all_dmrs <- function(dmr, arche, w, h) {
    dmr <- dmr[order(abs(dmr$meandiff), decreasing = TRUE),]
    dmr <- dmr[!is.na(dmr$overlapping.genes)]
    toPlot <- data.frame(
        gene = dmr$overlapping.genes,
        absmeandiff = abs(dmr$meandiff),
        fdr = dmr$HMFDR,
        rank = 1:length(dmr),
        direction = ifelse(dmr$meandiff>0, "Hypermethylation", "Hypomethylation")
    )
    # rename the long PCDHA family
    catch_me <- "PCDHA1, PCDHA2, PCDHA3, PCDHA4, PCDHA5, PCDHA6, PCDHA7, PCDHA8, PCDHA9, PCDHA10, PCDHA11, PCDHA12, PCDHA13, PCDHAC1"
    if (catch_me %in% toPlot$gene) toPlot$gene[toPlot$gene == catch_me] <- "PDCHA1-13, PCDHAC1"

    toPlot$rank <- factor(toPlot$rank, levels = toPlot$rank)

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/methylation/", arche, "_DMRs.png")
    png(filename, width = w, height = h, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = rank, y = absmeandiff, fill = direction)) +
        geom_col() +
        scale_x_discrete(labels = setNames(toPlot$gene, toPlot$rank)) +
        scale_fill_manual("Direction", values = binary_pal) +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),
            legend.position = c(0.99, 0.95),  
            legend.justification = c("right", "top")
        ) +
        labs(x = "Overlapping Genes", y = "Abs(Mean Methylation Difference)")
    )
    dev.off()
}

#' Plot DMRs using dmrcate
#' 
#' @param dmrs gr. output of extractRanges
#' @param annot output of cpg.annotate()
#' @param arche string
#' @param n int. DMR number in dmrs
#' @param gene string
#' @param cols string
#' @param h int height for figure
#' 
plot_dmr_arche <- function(dmrs, annot, arche, n, gene, h = 15) {
    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/methylation/DMR/", arche, "_", gene, ".png")
    png(filename, width = 18, height = h, res = 600, units = "in")
    print(
        DMR.plot(
        ranges=dmrs, dmr=n, CpGs=annot,
        what = "Beta", arraytype = "450K", phen.col=cols, genome = "hg19")
    )
    dev.off()
}

#' Plot overlap of genes in MYC gene sets and ARCHE2 DMRs
#' 
plot_arche2_myc_dmrs <- function(overlap, myc_targs) {

    # get MYC signatures with genes with DMRs
    myc <- data.frame(matrix(nrow=length(overlap), ncol=length(myc_targs)+1, data = 0))
    colnames(myc) <- c("Gene", names(myc_targs))
    myc$Gene <- overlap

    for (i in seq_along(overlap)) {
        gene <- overlap[i]
        for (sig in names(myc_targs)) {
            if (gene %in% myc_targs[[sig]]) myc[[sig]][i] <- 1
        }
    }
    myc <- myc[,-(which(colSums(myc[,-c(1)])==0)+1)] #offset by 1 to exclude "Gene"

    # get number of genes in each gene set
    gene_counts <- sapply(myc_targs, length)
    gene_counts <- data.frame(
        gene_set = names(myc_targs),
        count = gene_counts
    )

    # create toPlot
    toPlot <- reshape2::melt(myc, id.var = "Gene")
    toPlot$value <- factor(toPlot$value, levels = c(1, 0))
    gene_sets <- unique(toPlot$variable)
    toPlot$variable <- factor(toPlot$variable, levels = gene_sets)

    # heatmap
    p1 <- ggplot(toPlot, aes(x = Gene, y = variable, fill = value)) + 
        geom_tile(color = "gray") +
        scale_fill_manual(
            "Gene Present\nin Gene Set", 
            values = c(random_blue, "white"),
            labels = c("Present", "Absent")) +
        theme_void() +
        theme(
            axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),
            axis.text.y = element_text(size = 10, vjust = 0.5, hjust=1), 
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(y = "MYC Gene Set")

    # gene count
    toPlot <- gene_counts[rownames(gene_counts) %in% gene_sets,]
    toPlot$gene_set <- factor(toPlot$gene_set, levels = gene_sets)
    p2 <- ggplot() +
        geom_col(data = toPlot, aes(x = count, y = gene_set), fill = random_blue) +
        geom_text(data = toPlot, aes(x = count, y = gene_set, label = count), hjust = -0.1, size = 3) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
            ) +
        scale_x_continuous(limits = c(0, 1230), expand = c(0, 0)) +
        labs(y = NULL, x = "No. Genes\n\n")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    p1 <- p1+theme(legend.position = "none")

    filename <- "data/results/figures/2-MolecularSigAnalysis/methylation/ARCHE2_MYC_DMRs.png"
    png(filename, width = 9, height = 3, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, ncol = 3, nrow = 2,
        layout_matrix = rbind(c(1,1,2),
                            c(1,1,2)))
    )
    dev.off()

}