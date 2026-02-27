# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(survival)
    library(ggsurvfit)
    library(survminer)
})

source("utils/plots/signatures.R")
source("utils/get_data.R")
source("utils/palettes.R")

set.seed(123)

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)

# get signature scores from NMF
mat <- get_arche_tcga()
mat$variable <- meta$Sample.Name[match(mat$variable, meta$ATAC.Seq.File.Name)]

# read in tcga clinical data
pheno <- read.table("data/rawdata/tcga/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = TRUE, row.names = 1)
pheno <- t(pheno) |> as.data.frame()
rownames(pheno) <- gsub("\\.", "-", rownames(pheno))

###########################################################
# Get survival and PAM50 data
###########################################################

mat$OS <- pheno$overall_survival[match(mat$variable, rownames(pheno))] |> as.numeric()
mat$status <- pheno$status[match(mat$variable, rownames(pheno))] |> as.numeric()
mat$pam50 <- pheno$PAM50[match(mat$variable, rownames(pheno))] |> as.factor()

# clean dataframe
survdf <- unique(mat[,c("variable", "signature_assign", "subtype", "pam50", "status", "OS")])
colnames(survdf)[2] <- "ARCHE"

###########################################################
# Helper function to make survival plots
###########################################################

plot_survival <- function(df, ARCHE, label) {

    df$ARCHE <- ifelse(df$ARCHE == ARCHE, ARCHE, "Other")
    fit <- survfit(Surv(OS, status) ~ ARCHE, data = df)

    p <- ggsurvplot(
        fit,
        data = df,
        pval = TRUE,
        risk.table = TRUE,
        palette = c("#046C9A", "#C3BFCC"),
        legend.title = "ARCHE",
        xlab = "Time (days)"
    )
    p$plot <- p$plot + theme(
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)
    )

    filename <- paste0("data/results/figures/Misc/survivalplots/", label, ".png")
    png(filename, width=7, height=6, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()
}

###########################################################
# All samples and all ARCHEs
###########################################################

plot_survival(survdf, "ARCHE1", "all_ARCHE1")
plot_survival(survdf, "ARCHE2", "all_ARCHE2")
plot_survival(survdf, "ARCHE3", "all_ARCHE3")
plot_survival(survdf, "ARCHE4", "all_ARCHE4")
plot_survival(survdf, "ARCHE5", "all_ARCHE5")
plot_survival(survdf, "ARCHE6", "all_ARCHE6")

###########################################################
# Basal vs ARCHE2&5
###########################################################

basal <- survdf[survdf$subtype == "Basal",]

plot_survival(basal, "ARCHE5", "basal_ARCHE5")
plot_survival(basal, "ARCHE2", "basal_ARCHE2")

###########################################################
# Luminals vs ARCHE1,4,&6 together
###########################################################

lum <- survdf[survdf$subtype %in% c("LumA", "LumB"),]

plot_survival(lum, "ARCHE1", "lum_ARCHE1")
plot_survival(lum, "ARCHE4", "lum_ARCHE4")
plot_survival(lum, "ARCHE6", "lum_ARCHE6")

###########################################################
# LuminalA vs ARCHE1,4,&6
###########################################################

lumA <- survdf[survdf$subtype == "LumA",]

plot_survival(lumA, "ARCHE1", "lumA_ARCHE1")
plot_survival(lumA, "ARCHE4", "lumA_ARCHE4")
plot_survival(lumA, "ARCHE6", "lumA_ARCHE6")

###########################################################
# LuminalB vs ARCHE1,4,&6
###########################################################

lumB <- survdf[survdf$subtype == "LumB",]

plot_survival(lumB, "ARCHE1", "lumB_ARCHE1") # !!!!!!
plot_survival(lumB, "ARCHE4", "lumB_ARCHE4")
plot_survival(lumB, "ARCHE6", "lumB_ARCHE6")

###########################################################
# LuminalB vs ARCHE1 and all others
###########################################################

lumB <- survdf[survdf$subtype == "LumB",]

fit <- survfit(Surv(OS, status) ~ ARCHE, data = lumB)

p <- ggsurvplot(
    fit,
    data = lumB,
    pval = TRUE,
    risk.table = TRUE,
    palette = unname(ARCHE_pal)[c(1:4,6)],
    legend.title = "ARCHE",
    xlab = "Time (days)"
)
p$plot <- p$plot + theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
)

filename <- paste0("data/results/figures/Misc/survivalplots/", "lumb_all_ARCHEs", ".png")
png(filename, width=8, height=8, units='in', res = 600, pointsize=80)
print(p)
dev.off()

###########################################################
# All ARCHEs together 
###########################################################

fit <- survfit(Surv(OS, status) ~ ARCHE, data = survdf)

p <- ggsurvplot(
    fit,
    data = survdf,
    pval = TRUE,
    risk.table = TRUE,
    palette = unname(ARCHE_pal),
    legend.title = "ARCHE",
    xlab = "Time (days)"
)
p$plot <- p$plot + theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
)

filename <- paste0("data/results/figures/Misc/survivalplots/", "all_ARCHEs", ".png")
png(filename, width=8, height=8, units='in', res = 600, pointsize=80)
print(p)
dev.off()