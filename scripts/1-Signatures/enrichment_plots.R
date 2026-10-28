# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(viridis)
    library(patchwork)
    library(rrvgo)
})

source("utils/plots/signatures.R")
source("utils/palettes.R")

set.seed(123)

###########################################################
# Load in GREAT enrichment
###########################################################

INDIR <- "data/results/data/1-Signatures/GREAT/"

great1 <- fread(paste0(INDIR, "ARCHE1_50k.tsv"), data.table = FALSE)
great2 <- fread(paste0(INDIR, "ARCHE2_50k.tsv"), data.table = FALSE)
great3 <- fread(paste0(INDIR, "ARCHE3_50k.tsv"), data.table = FALSE)
great4 <- fread(paste0(INDIR, "ARCHE4_50k.tsv"), data.table = FALSE)
great5 <- fread(paste0(INDIR, "ARCHE5_50k.tsv"), data.table = FALSE)
great6 <- fread(paste0(INDIR, "ARCHE6_50k.tsv"), data.table = FALSE)

###########################################################
# Format for plotting
###########################################################

# helper function
format_great <- function(df, arche) {
    df <- df[df$Label == "Biological Process",] # keep only biological process
    df$ARCHE <- arche
    df$pair <- paste0(df$name, "_", df$ARCHE)
    df <- df[df$Total_Genes_Annotated > 5,]
    return(df)
}

great1 <- format_great(great1, "ARCHE1")
great2 <- format_great(great2, "ARCHE2")
great3 <- format_great(great3, "ARCHE3")
great4 <- format_great(great4, "ARCHE4")
great5 <- format_great(great5, "ARCHE5")
great6 <- format_great(great6, "ARCHE6")


top_pairs <- c(
    great1$pair[1:10],
    great2$pair[1:10],
    great3$pair[1:10],
    great4$pair[1:10],
    great5$pair[1:10],
    great6$pair[1:10]
)

toPlot <- rbind(great1, great2, great3, great4, great5, great6)
toPlot <- toPlot[toPlot$name %in% sub("_.*", "", top_pairs),]

# shorten names
toPlot$name <- sub(
    "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",
    "oxidoreductase activity, acting on CH-OH donor group, NAD(P) as acceptor",
    toPlot$name
)
toPlot$name <- sub(
    "oxidoreductase activity, acting on CH-OH group of donors",
    "oxidoreductase activity, acting on CH-OH donor group",
    toPlot$name
)

toPlot$ARCHE <- factor(toPlot$ARCHE, levels = paste0("ARCHE", 6:1))

toPlot <- toPlot[complete.cases(toPlot$name),]

###########################################################
# Read in annotation
###########################################################

anno <- read.csv("utils/anno/GO_term_categories.csv")
toPlot$Category <- anno$category[match(toPlot$name, anno$name)]

anno <- anno[order(anno$category),]
toPlot$name <- factor(toPlot$name, levels = unique(anno$name))

###########################################################
# Plot top GO terms
###########################################################

p <- ggplot(toPlot, aes(x = name, y = ARCHE, size = -log10(Hyper_Adjp_BH + 0.001), color = Hyper_Fold_Enrichment)) +
    annotate("rect", ymin = -Inf, ymax = Inf, 
        xmin = c(0.5, 7.5, 23.5, 27.5, 40.5, 48.5, 54.5),
        xmax = c(4.5, 16.5, 26.5, 28.5, 45.5, 50.5, 55.5),
        fill = "#CCCFD5", alpha = 0.5) +
    geom_point() +
    scale_color_viridis_c(option = "mako", direction = -1, end = 0.9) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),
        axis.title.y = element_blank()
    ) +
    labs(size = "-log(FDR)", color = "Fold\nEnrichment", x = "GO Term")


ggsave("data/results/figures/1-Signatures/GREAT/GREAT_50k_top10BP.png", p, width = 10.5, height = 6)

###########################################################
# Collapse by rrvgo
###########################################################

# filter significant BP terms (10-500 annotated genes and fold enrichment >= 1.3)
great_all <- rbind(great1, great2, great3, great4, great5, great6)
great_all <- great_all[
    great_all$Total_Genes_Annotated >= 10 &
    great_all$Total_Genes_Annotated <= 500 &
    great_all$Hyper_Fold_Enrichment >= 1.3,
]

# calculate similarity across terms
sim <- calculateSimMatrix(
    unique(great_all$ID),
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel"
)

# score by best FDR across ARCHEs and reduce redundant terms
scores <- tapply(-log10(great_all$Hyper_Adjp_BH), great_all$ID, max)[rownames(sim)]
reduced <- reduceSimMatrix(sim, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")

# assign terms to rrvgo cluster
great_all$cluster <- reduced$cluster[match(great_all$ID, reduced$go)]
great_all$parentTerm <- reduced$parentTerm[match(great_all$ID, reduced$go)]
great_all <- great_all[!is.na(great_all$parentTerm),] # terms without IC information in GO

# summarise clusters
clusters <- merge(
    aggregate(ID ~ cluster + parentTerm + ARCHE, data = great_all, length),
    aggregate(cbind(Hyper_Fold_Enrichment, Hyper_Adjp_BH) ~ cluster + parentTerm + ARCHE, data = great_all,
        FUN = function(x) c(mean = mean(x), min = min(x))),
    by = c("cluster", "parentTerm", "ARCHE")
)
clusters <- data.frame(
    cluster = clusters$cluster,
    parentTerm = clusters$parentTerm,
    ARCHE = clusters$ARCHE,
    n_terms = clusters$ID,
    mean_FE = clusters$Hyper_Fold_Enrichment[, "mean"],
    best_FDR = clusters$Hyper_Adjp_BH[, "min"]
)

# count number of ARCHEs for each cluster
clusters$n_ARCHEs <- ave(clusters$cluster, clusters$cluster, FUN = length)
clusters <- clusters[order(clusters$ARCHE, -clusters$n_terms, clusters$best_FDR),]

write.csv(great_all, "data/results/data/1-Signatures/GREAT/rrvgo_terms.csv", row.names = FALSE)
write.csv(clusters, "data/results/data/1-Signatures/GREAT/rrvgo_clusters.csv", row.names = FALSE)

###########################################################
# Select ARCHE-specific rrvgo clusters
###########################################################

# keep clusters enriched in <=3 ARCHEs
specific <- clusters[clusters$n_ARCHEs <= 3,]

# assign cluster to its most significant ARCHE
owner <- specific[order(specific$best_FDR),]
owner <- owner[!duplicated(owner$cluster), c("cluster", "ARCHE")]
colnames(owner)[2] <- "owner"

# top 3 clusters per ARCHE by FDR
top_clusters <- merge(specific, owner, by = "cluster")
top_clusters <- top_clusters[top_clusters$ARCHE == top_clusters$owner,]
top_clusters <- top_clusters[order(top_clusters$ARCHE, top_clusters$best_FDR),]
top_clusters <- do.call(rbind, lapply(split(top_clusters, top_clusters$ARCHE), head, 3))

# rename cluster by most significant term
owned_terms <- great_all[paste(great_all$cluster, great_all$ARCHE) %in% paste(top_clusters$cluster, top_clusters$ARCHE),]
owned_terms <- owned_terms[order(owned_terms$Hyper_Adjp_BH),]
owned_terms <- owned_terms[!duplicated(owned_terms$cluster),]
top_clusters$label <- owned_terms$name[match(top_clusters$cluster, owned_terms$cluster)]

# shorten or clarify labels (using GO terms within the same cluster)
relabel <- c(
    "negative regulation of pri-miRNA transcription from RNA polymerase II promoter" = "negative regulation of miRNA transcription",
    "positive regulation of metanephros development" = "regulation of metanephros development",
    "regulation of cell activation" = "regulation of leukocyte activation"
)
top_clusters$label <- ifelse(top_clusters$label %in% names(relabel), relabel[top_clusters$label], top_clusters$label)

###########################################################
# Plot ARCHE-specific rrvgo clusters
###########################################################

# all ARCHEs enriched for the selected clusters
toPlot <- specific[specific$cluster %in% top_clusters$cluster,]
toPlot$label <- top_clusters$label[match(toPlot$cluster, top_clusters$cluster)]

# order clusters along the diagonal (ARCHE1 -> ARCHE6) and facet by owning ARCHE
toPlot$label <- factor(toPlot$label, levels = rev(top_clusters$label))
toPlot$ARCHE <- factor(toPlot$ARCHE, levels = paste0("ARCHE", 1:6))
toPlot$group <- factor(top_clusters$ARCHE[match(toPlot$cluster, top_clusters$cluster)], levels = paste0("ARCHE", 1:6))
toPlot$group <- droplevels(toPlot$group)

p <- ggplot(toPlot, aes(x = ARCHE, y = label, color = mean_FE, size = -log10(best_FDR))) +
    geom_point() +
    ggh4x::facet_grid2(
        group ~ ., scales = "free_y", space = "free_y",
        strip = ggh4x::strip_themed(
            background_y = ggh4x::elem_list_rect(fill = ARCHE_pal[levels(toPlot$group)])
        )
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_color_viridis_c("Mean Fold\nEnrichment", option = "mako", direction = -1, end = 0.9) +
    scale_size("-log10(FDR)", range = c(1.5, 6)) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(size = 0, margin = margin(0, 1, 0, 1, "mm")), # narrow coloured strip, no label
        panel.spacing = unit(1, "mm"),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 8),
        panel.grid.minor = element_blank()
    )

ggsave("data/results/figures/1-Signatures/GREAT/GREAT_50k_rrvgo_collapsed.png", p, width = 5, height = 4)
