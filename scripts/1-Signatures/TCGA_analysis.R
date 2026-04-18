# demonstrate BCa representativeness from TCGA

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
})

source("utils/palettes.R")

set.seed(123)

OUTDIR <- "data/results/figures/1-Signatures/TCGA_ARCHE_vars/"

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)

# read in tcga clinical data
pheno <- read.table("data/rawdata/tcga/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = TRUE, row.names = 1)
pheno <- t(pheno) |> as.data.frame()
rownames(pheno) <- gsub("\\.", "-", rownames(pheno))

###########################################################
# Merge clinical variables and format
###########################################################

for (column in colnames(pheno)) {
    meta[[column]] <- pheno[[column]][match(meta$Sample.Name, rownames(pheno))]
}

# format clinical vars
meta$Subtype <- factor(meta$Subtype, levels = rev(names(subtype_pal)))
meta$Tumor_purity <- as.numeric(meta$Tumor_purity)
meta$years_to_birth <- as.numeric(meta$years_to_birth)
meta$overall_survival <- as.numeric(meta$overall_survival)
meta$pathologic_stage <- factor(meta$pathologic_stage, levels = names(stage_pal))

###########################################################
# Plot subtype
###########################################################

p <- ggplot(meta, aes(y = Subtype, fill = Subtype)) +
  geom_bar(color = "black", width = 0.95, linewidth = 0.3) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count), colour = Subtype %in% c("Basal", "LumB")),
    x = 1.25,
    hjust = 0,
    vjust = 0.5,
    size = 3
  ) +
  scale_fill_manual(values = subtype_pal) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10)) +
  labs(x = "Tumour Sample Count")

ggsave(paste0(OUTDIR, "subtype.png"), p, width = 3, height = 2.5, bg = "transparent")

###########################################################
# Plot stage and ARCHE score
###########################################################

# helper function to plot proportions
plot_prop <- function(column_name, legend_title, pal, filename) {
    p <- ggplot(meta, aes(y = Subtype, fill = .data[[column_name]])) +
        geom_bar(position = "fill", color = "black", width = 0.95, linewidth = 0.3) +
        scale_x_continuous(labels = scales::percent, limits = c(-0.05, 1.05)) +
        scale_fill_manual(values = pal) +
        theme_bw() +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10)
        ) +
        labs(x = "Proportion", fill = legend_title)

    ggsave(paste0(OUTDIR, filename), p, width = 3, height = 2.5, bg = "transparent")
}

plot_prop("pathologic_stage", "Stage", stage_pal, "stage.png")
plot_prop("ARCHE", "ARCHE", ARCHE_pal, "ARCHE.png")

###########################################################
# Plot distribution of other clinical variables
###########################################################

# helper function to plot boxplots
plot_boxplot <- function(column_name, axis_name, filename) {

    p <- ggplot(meta, aes(y = Subtype, x = .data[[column_name]], fill = Subtype)) +
        geom_boxplot() +
        geom_jitter(height = 0.1) +
        scale_fill_manual(values = subtype_pal) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10)
        ) +
        labs(x = axis_name)

    ggsave(paste0(OUTDIR, filename), p, width = 2, height = 2.5, bg = "transparent")

}

plot_boxplot("Tumor_purity", "Tumour Purity", "purity.png")
plot_boxplot("years_to_birth", "Age", "age.png")
plot_boxplot("overall_survival", "OS (Days)", "os.png")