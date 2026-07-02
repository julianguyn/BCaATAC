# demonstrate BCa representativeness from TCGA

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(survival)
    library(ggsurvfit)
    library(survminer)
    library(patchwork)
})

source("utils/palettes.R")

set.seed(123)

# set variables based on analysis
args <- commandArgs(trailingOnly = TRUE)
ref <- args[1] # either "Subtype" or "ARCHE"

OUTDIR <- paste0("data/results/figures/1-Signatures/TCGA_ARCHE_vars/", ref, "_")
main_pal <- switch(ref, Subtype = subtype_pal, ARCHE = ARCHE_pal)
comp_pal <- switch(ref, Subtype = ARCHE_pal, ARCHE = subtype_pal)

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

# load in TMB
tmb <- readRDS("data/procdata/TCGA/tmb.rds")

###########################################################
# Merge clinical variables and format
###########################################################

for (column in colnames(pheno)) {
    meta[[column]] <- pheno[[column]][match(meta$Sample.Name, rownames(pheno))]
}

# format clinical vars
meta$Subtype <- factor(meta$Subtype, levels = rev(names(subtype_pal)))
meta$ARCHE <- factor(meta$ARCHE, levels = rev(names(ARCHE_pal)))
meta$Tumor_purity <- as.numeric(meta$Tumor_purity)
meta$years_to_birth <- as.numeric(meta$years_to_birth)
meta$overall_survival <- as.numeric(meta$overall_survival)
meta$status <- as.numeric(meta$status)
meta$pathologic_stage <- factor(meta$pathologic_stage, levels = rev(names(stage_pal)))

###########################################################
# Plot subtype/ARCHE
###########################################################

p <- ggplot(meta, aes(y = .data[[ref]], fill = .data[[ref]])) +
  geom_bar(color = "black", width = 0.95, linewidth = 0.3) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count), colour = .data[[ref]] %in% c("Basal", "LumB")),
    x = 1.25,
    hjust = 0,
    vjust = 0.5,
    size = 3
  ) +
  scale_fill_manual(values = main_pal) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10)) +
  labs(x = "Tumour Sample Count")

ggsave(paste0(OUTDIR, ".png"), p, width = 3, height = 2.5, bg = "transparent")

###########################################################
# Plot stage and ARCHE score
###########################################################

# helper function to plot proportions
plot_prop <- function(column_name, legend_title, pal, filename) {
    p <- ggplot(meta, aes(y = .data[[ref]], fill = .data[[column_name]])) +
        geom_bar(position = "fill", color = "black", width = 0.95, linewidth = 0.3) +
        scale_x_continuous(labels = scales::percent, limits = c(-0.05, 1.05)) +
        scale_fill_manual(values = pal) +
        theme_bw() +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10),
            legend.position = "none"
        ) +
        labs(x = "Proportion", fill = legend_title)

    ggsave(paste0(OUTDIR, filename), p, width = 2, height = 2.5, bg = "transparent")
}

plot_prop("pathologic_stage", "Stage", stage_pal, "stage.png")
if (ref == "Subtype") {
    plot_prop("ARCHE", "ARCHE", comp_pal, "ARCHE.png")
} else {
    plot_prop("Subtype", "Subtype", comp_pal, "subtype.png")
}

###########################################################
# Plot distribution of other clinical variables
###########################################################

# helper function to plot boxplots
plot_boxplot <- function(column_name, axis_name, filename) {

    p <- ggplot(meta, aes(y = .data[[ref]], x = .data[[column_name]], fill = .data[[ref]])) +
        geom_boxplot() +
        geom_jitter(height = 0.1) +
        scale_fill_manual(values = main_pal) +
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

###########################################################
# Plot TMB box plots
###########################################################

tmb$Subtype <- factor(tmb$Subtype, levels = rev(names(subtype_pal)))
tmb$ARCHE <- factor(tmb$ARCHE, levels = rev(names(ARCHE_pal)))

p <- ggplot(tmb, aes(y = .data[[ref]], x = total_perMB_log, fill = .data[[ref]])) +
        geom_boxplot() +
        geom_jitter(height = 0.1) +
        scale_fill_manual(values = main_pal) +
        scale_x_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, by = 0.4)) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 10)
        ) +
        labs(x = "TMB/MB (log10)")
ggsave(paste0(OUTDIR, "tmb.png"), p, width = 2, height = 2.5, bg = "transparent")

###########################################################
# Plot survival plots
###########################################################

# cutoff times
cutoff <- 1095
meta$time_cutoff <- pmin(meta$overall_survival, cutoff)
meta$event_cutoff <- ifelse(meta$overall_survival <= cutoff, meta$status, 0)

# plot survival plots
plot_survival <- function(df, cutoff) {

    df$ref <- ifelse(df[[ref]] == group, group, "Other")
    df$ref <- factor(df$ref, levels = c("Other", group))
    fit <- survfit(Surv(time_cutoff, event_cutoff) ~ ref, data = df)

    p <- ggsurvplot(
        fit,
        data = df,
        size = 2,
        pval = TRUE,
        risk.table = TRUE,
        palette = c("#979797", main_pal[[group]]),
        legend.title = "",
        xlab = "Time (days)",
        xlim = c(0, cutoff),
        break.time.by = 365
    )
    p$plot <- p$plot + 
        theme(
            legend.key.size = unit(3, "cm"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            
        )
    p$table <- p$table + theme(plot.title = element_text(size = 12))

    filename <- paste0("data/results/figures/1-Signatures/survivalplots/", ref, "_", group, ".png")
    png(filename, width=4.5, height=5, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()
}


for (group in unique(df[[ref]])) {
    plot_survival(meta, cutoff)
}




# --- sanity check of survival by subtypes across ALL BRCA tumours
pheno <- pheno[-which(is.na(pheno$PAM50)),] #826
pheno$overall_survival <- as.numeric(pheno$overall_survival)
pheno$status <- as.numeric(pheno$status)
pheno$PAM50 <- as.character(pheno$PAM50)

for (group in unique(pheno$PAM50)) {
    pheno$ref <- ifelse(pheno$PAM50 == group, group, "Other")
    fit <- survfit(Surv(overall_survival, status) ~ ref, data = pheno)

    p <- ggsurvplot(
        fit,
        data = pheno,
        pval = TRUE,
        risk.table = TRUE,
        palette = c("#046C9A", "#827A6F"),
        legend.title = "",
        xlab = "Time (days)",
        xlim = c(0, 1825),
        break.time.by = 365
    )
    p$plot <- p$plot + theme(
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)
    )

    filename <- paste0("data/results/figures/1-Signatures/survivalplots/ALL_BRCA_", group, ".png")
    png(filename, width=7, height=6, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()
}

fit <- survfit(Surv(overall_survival, status) ~ PAM50, data = pheno)

p <- ggsurvplot(
    fit,
    data = pheno,
    pval = TRUE,
    risk.table = TRUE,
    palette = unname(subtype_pal[1:4]),
    legend.title = "",
    xlab = "Time (days)",
    xlim = c(0, 1825),
    break.time.by = 365
)
p$plot <- p$plot + theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
)

filename <- paste0("data/results/figures/1-Signatures/survivalplots/ALL_BRCA_all_subtypes.png")
png(filename, width=7, height=6, units='in', res = 600, pointsize=80)
print(p)
dev.off()