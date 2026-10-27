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

# load in pam50 subtyping
pam50 <- readRDS("data/procdata/TCGA/pam50_subtyping_full_cohort.rds")
pam50_scores <- as.data.frame(pam50$subtype.proba)
pam50_scores$assigned <- as.character(pam50$subtype)
rownames(pam50_scores) <- gsub("\\.", "-", rownames(pam50_scores))

###########################################################
# Merge clinical variables and format
###########################################################

for (column in colnames(pheno)) {
    meta[[column]] <- pheno[[column]][match(meta$Sample.Name, rownames(pheno))]
}
meta$Subtype <- pam50_scores$assigned[match(meta$Sample.Name, rownames(pam50_scores))]

# format clinical vars
meta$Subtype <- factor(meta$Subtype, levels = names(subtype_pal))
meta$ARCHE <- factor(meta$ARCHE, levels = names(ARCHE_pal))
meta$Tumor_purity <- as.numeric(meta$Tumor_purity)
meta$years_to_birth <- as.numeric(meta$years_to_birth)
meta$overall_survival <- as.numeric(meta$overall_survival)
meta$status <- as.numeric(meta$status)
meta$pathologic_stage <- factor(meta$pathologic_stage, levels = rev(names(stage_pal)))

###########################################################
# Plot subtype/ARCHE
###########################################################

p4 <- ggplot(meta, aes(x = .data[[ref]], fill = .data[[ref]])) +
  geom_bar(color = "black", width = 0.95, linewidth = 0.3) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count), colour = .data[[ref]] %in% c("Basal", "LumB", "Normal", "ARCHE1", "ARCHE5")),
    y = 2.5,
    size = 3
  ) +
  scale_y_continuous(limits = c(0, 24), expand = c(0,0)) +
  scale_fill_manual(values = main_pal) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
    axis.title.y = element_text(size = 10)) +
  labs(y = "No. Samples")

###########################################################
# Plot distribution of other clinical variables
###########################################################

# helper function to plot boxplots
plot_boxplot <- function(column_name, axis_name, filename) {

    p <- ggplot(meta, aes(x = .data[[ref]], y = .data[[column_name]], fill = .data[[ref]])) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.6) +
        scale_fill_manual(values = main_pal) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 9)
        ) +
        labs(y = axis_name) 
    return(p)
}

p1 <- plot_boxplot("Tumor_purity", "Tumour Purity", "purity.png") + ylim(c(0.4, 1.05))
p2 <- plot_boxplot("years_to_birth", "Age", "age.png")

###########################################################
# Plot TMB box plots
###########################################################

tmb$Subtype <- pam50_scores$assigned[match(tmb$Sample_ID, rownames(pam50_scores))]
tmb$Subtype <- factor(tmb$Subtype, levels = names(subtype_pal))
tmb$ARCHE <- factor(tmb$ARCHE, levels = names(ARCHE_pal))

p3 <- ggplot(tmb, aes(y = .data[[ref]], x = total_perMB_log, fill = .data[[ref]])) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(height = 0.1, alpha = 0.6) +
        scale_fill_manual(values = main_pal) +
        scale_x_continuous(limits = c(-0.85, 1), breaks = seq(-0.8, 0.8, by = 0.4)) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 9)
        ) +
        labs(x = "TMB/MB (log10)") +
        coord_flip()


###########################################################
# ANOVA and Tukey
###########################################################

TukeyHSD(aov(Tumor_purity ~ ARCHE, data = meta)) # largest ARCHE5 p is 0.006
TukeyHSD(aov(years_to_birth ~ ARCHE, data = meta))
TukeyHSD(aov(total_perMB_log ~ ARCHE, data = tmb)) # ARCHE3-6 p is 0.0097

TukeyHSD(aov(Tumor_purity ~ Subtype, data = meta))
TukeyHSD(aov(years_to_birth ~ Subtype, data = meta))
TukeyHSD(aov(total_perMB_log ~ Subtype, data = tmb)) # LumA-Basal p = 0.021, LumA-HER2 p = 0.00004, LumB-HER2 p = 0.04 

if (ref == "ARCHE") {
    p1 <- p1 + annotate("text", x = "ARCHE5", y = 1, label = "**", size = 4)
    p3 <- p3 + geom_signif(y_position = 0.77, xmin = 3, xmax = 6, annotation = "*", tip_length = 0, textsize = 4, vjust = 0.55)
} else {
    p3 <- p3 +
        geom_signif(y_position = 0.7, xmin = 2, xmax = 3, annotation = "**", tip_length = 0, textsize = 4, vjust = 0.55) +
        geom_signif(y_position = 0.84, xmin = 1, xmax = 3, annotation = "*", tip_length = 0, textsize = 4, vjust = 0.55) +
        geom_signif(y_position = 0.98, xmin = 2, xmax = 4, annotation = "*", tip_length = 0, textsize = 4, vjust = 0.55)
}

###########################################################
# Merge all plots
###########################################################

p <- p1 / p3 / p2 / p4 + plot_layout(height = c(2,2,2,1))
w <- switch(ref, ARCHE = 2, Subtype = 1.8)
h <- switch(ref, ARCHE = 6, Subtype = 5.9)
ggsave(paste0(OUTDIR, "boxplots.png"), p, width = w, height = h, bg = "transparent")

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
