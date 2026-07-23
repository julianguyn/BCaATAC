# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(viridis)
    library(patchwork)
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
# Collapse by category
###########################################################

df <- aggregate(Hyper_Fold_Enrichment ~ Category + ARCHE, data = toPlot, sum)
counts <- aggregate(ID ~ Category + ARCHE, data = toPlot, length)
df$count <- counts$ID
df$ARCHE <- factor(df$ARCHE, levels = paste0("ARCHE", 1:6))

p <- ggplot(df, aes(x = ARCHE, y = Category, color = Hyper_Fold_Enrichment, size = count)) +
    geom_point() +
    scale_color_viridis_c("Sum of\nHyper Fold\nEnrichment", option = "mako", direction = -1, end = 0.9) +
    scale_size("Count of\nGO Terms", range = c(2, 10)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        legend.key.size = unit(0.4, "cm")
    )
ggsave("data/results/figures/1-Signatures/GREAT/GREAT_50k_top10BP_collapsed.png", p, width = 5, height = 4.5)
