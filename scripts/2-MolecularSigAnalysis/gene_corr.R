# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(viridis)
})

source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load in counts matrix
counts <- get_tcga_rna()
counts <- round(2^counts - 1) # unlog2normalize

# load in metadata file
meta <- read.csv("data/rawdata/tcga/TCGA_sourcefiles.csv")
meta <- meta[match(colnames(counts), meta$Sample.Name),]

# get signature scores from NMF
mat <- get_arche_tcga()
mat$variable <- meta$Sample.Name[match(mat$variable, meta$ATAC.Seq.File.Name)]

# format into wide matrix
arche_df <- mat[,1:3] %>%
    pivot_wider(
        names_from = variable,
        values_from = value
    ) %>% as.data.frame()
rownames(arche_df) <- arche_df$ARCHE
arche_df$ARCHE <- NULL
arche_df <- arche_df[,match(colnames(counts), colnames(arche_df))]

###########################################################
# Gene expression correlation
###########################################################

genes  <- rownames(counts)
arches <- rownames(arche_df)

results_df <- data.frame(
    ARCHE = rep(arches, each = length(genes)),
    gene = rep(genes, length(arches)),
    rho = NA,
    pvalue = NA
)

counter <- 1
for (gene in genes) {
    for (arche in arches) {
        
        ct <- cor.test(
        as.numeric(counts[gene, ]),
        as.numeric(arche_df[arche, ]),
        method = "spearman"
        )

        results_df$rho[counter] <- ct$estimate
        results_df$pvalue[counter] <- ct$p.value
        
        counter <- counter + 1
    }
}

results_df$FDR <- p.adjust(results_df$pvalue, method = "BH", n = length(results_df$pvalue))
saveRDS(results_df, file = "data/results/data/2-MolecularSigAnalysis/gene_correlations_arches.rds")

###########################################################
# Identify top genes and plot
###########################################################

sig_res <- results_df[results_df$FDR < 0.05,]
sig_res <- sig_res[order(abs(sig_res$rho), decreasing = TRUE),]

# plot top 10 per ARCHE
for (arche in paste0("ARCHE", 1:6)) {
    toPlot <- sig_res[sig_res$ARCHE == arche,][1:10,]
    toPlot$gene <- factor(toPlot$gene, levels = rev(toPlot$gene))
    
    p <- ggplot(toPlot, aes(x = "Spearman rho", y = gene, fill = rho, size = -log2(FDR))) +
        geom_point(shape = 21) +
        geom_text(aes(label = round(rho, 2)), size = 2) +
        scale_fill_gradient2(high = binary_pal[1], low = binary_pal[2], midpoint = 0) +
        scale_size_continuous(range = c(5, 10)) +
        theme_bw() +
        theme(
            axis.title.y = element_blank(), axis.title.x = element_blank()
        ) +
        ggtitle(arche)

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/gene_corr/", arche, ".png")
    ggsave(filename, p, width = 3, height = 4)
}
