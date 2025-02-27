setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(meta)
})

###########################################################
# Load in data
###########################################################

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# load in drug response associations
load("DrugResponse/results/data/indiv_PSet_CI.RData")
all_com <- rbind(ubr1_com, ubr2_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com) #combine all drug response from all psets


###########################################################
# Identify Class B Biomarkers
###########################################################

# keep only signature-drug pairs that are significant in 2+ PSets
df <- all_com[which(all_com$pairs %in% unique(sig_com$pairs[duplicated(sig_com$pairs)])),]
df$sig <- ifelse(df$ci > 0.6 | df$ci < 0.4, "sig", "not sig")
df[df$FDRsig == FALSE,]$sig <- "not sig"
df$pset <- toupper(df$pset)

# specify sensitivity and resistance
df$type <- ifelse(df$ci > 0.5, "Sensitivity", "Resistance")
df$type <- factor(df$type, levels = c("Sensitivity", "Resistance"))

# identify class B
df$classB <- FALSE
for (pair in unique(df$pairs)) {
    # for each pair, check if sig associations are all the same type (all sensitivity or all resistance)
    tmp <- df[df$pairs == pair,]
    ifelse(length(unique(tmp[tmp$sig == "sig",]$type)) == 1, df[df$pairs == pair,]$classB <- TRUE, df[df$pairs == pair,]$classB <- FALSE)
}

###########################################################
# Define palettes for plotting
###########################################################

# signature palette
pal <- c("CAS-1" = "#046C9A", "CAS-2" = "#BBADB9", "CAS-3" = "#7294D4", 
        "CAS-4" = "#E8E1D9", "CAS-5" = "#AFC5D8", "CAS-6" = "#DF9C93")
# sensitive vs resistant palette
pal2 <- c("#899DA4", "#BC4749")

###########################################################
# Plot all biomarkers in >=2 PSets
###########################################################

png("DrugResponse/results/figures/robust_PSet/two_psets_all.png", width = 17, height = 5, res = 600, units = "in")
ggplot(df, aes(x = pset, y = ci - 0.5, fill = ifelse(sig == "sig", type, "Not Significant"))) + geom_bar(stat="identity", color = "black") +
    facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0), labels = function(y) y + 0.5, oob = scales::squish) +
    scale_fill_manual(values = c(pal2, "#E8E1D9"), labels = c("Sensitivity", "Resistance", "Not Significant")) +
    labs(fill = "CAS", y = "Concordance Index (CI)", x = "PSet") + 
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") + geom_hline(yintercept = c(0)) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


###########################################################
# Save Class B Biomarkers
###########################################################

# remove non-significant associations
df <- df[df$classB == TRUE,]
df <- df[df$sig == "sig",]

write.csv(df, file = "DrugResponse/results/data/ClassB_Biomarkers.csv", quote = F, row.names = F)

###########################################################
# Plot Class B Biomarkers
###########################################################

png("DrugResponse/results/figures/robust_PSet/classB.png", width = 14, height = 5, res = 600, units = "in")
ggplot(df, aes(x = pset, y = ci - 0.5, fill = signature)) + geom_bar(stat="identity", color = "black") +
    facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0), labels = function(y) y + 0.5, oob = scales::squish) +
    scale_fill_manual(values = pal) +
    labs(fill = "CAS", y = "Concordance Index (CI)", x = "PSet") + 
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") + geom_hline(yintercept = c(0)) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
          axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 16), 
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 16),
          strip.text = element_text(size = 13),
          axis.title.y = element_text(size = 16))
dev.off()


###########################################################
# Perform Meta-Analysis
###########################################################

# keep only signature-drug pairs that are present in at least 3 PSets
df <- all_com[which(all_com$pairs %in% names(table(all_com$pairs)[table(all_com$pairs) > 2])),]
df$sig <- ifelse(df$ci > 0.6 | df$ci < 0.4, "sig", "not sig")
df[df$FDRsig == FALSE,]$sig <- "not sig"
df$pset <- toupper(df$pset)

# data frame to hold meta estimates
estimates <- as.data.frame(matrix(data = NA, nrow = length(unique(df$pairs)), ncol = 9))
colnames(estimates) <- c("signature","drug", "TE", "seTE", "upper", "lower", "pval", "pval.Q", "I2")

# perform meta-analysis
for (i in 1:length(unique(df$pair))) {
  print(paste(i, "out of", length(unique(df$pair))))
  
  pair <- unique(df$pairs)[i]
  tmp <- df[which(df$pair == pair),]
  meta <- metagen(ci, se, data = tmp, method.tau = "DL", studlab = tmp$pset)
  
  estimates$signature[i] <- tmp$signature[1]
  estimates$drug[i] <- tmp$drug[1]
  estimates$TE[i] <- meta$TE.random
  estimates$seTE[i] <- meta$seTE.random
  estimates$upper[i] <- meta$upper.random
  estimates$lower[i] <- meta$lower.random
  estimates$pval[i] <- meta$pval.random
  estimates$pval.Q[i] <- meta$pval.Q
  estimates$I2[i] <- meta$I2
  
  # Forest Plots
  fileName = paste0("DrugResponse/results/figures/robust_PSet/meta/",df$signature[i],"_",df$drug[i],".png")
  
  png(fileName, width = 10, height = 4, res = 600, units = "in")
  title <- paste0(df$circ_ids[i], "/", df$drug[i])
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

estimates$FDR <- p.adjust(estimates$pval, method = "BH", n = length(estimates$pval))
estimates$pairs <- paste0(estimates$signature, "_", estimates$drug)

###########################################################
# Save Meta-Analysis estimates
###########################################################

write.csv(estimates, file = "DrugResponse/results/tables/robust_PSet/meta_estimates.csv", quote = F, row.names = F)


###########################################################
# Identify Class C Biomarkers
###########################################################

classC <- estimates[which(estimates$TE > 0.6 | estimates$TE < 0.4),]
classC <- classC[classC$FDR < 0.05,]
classC$type <- ifelse(classC$TE > 0.5, "Sensitivy", "Resistance")

###########################################################
# Save Class C Biomarkers
###########################################################

write.csv(classC, file = "DrugResponse/results/data/ClassC_Biomarkers.csv", quote = F, row.names = F)

###########################################################
# Format dataframe to plot Class C biomarkers
###########################################################

# only plot classC
df <- df[which(df$pairs %in% classC$pairs),]
df$FDRsig <- factor(df$FDRsig, levels = c(TRUE, FALSE))
df$meta <- factor(ifelse(df$pset == "Meta Estimate", TRUE, FALSE), levels = c(TRUE, FALSE))

# add meta estimate to df
df <- rbind(df, data.frame(signature = classC$signature, drug = classC$signature, 
                            ci = classC$TE, pvalue = classC$pval, se = classC$seTE,
                            upper = classC$upper, lower = classC$lower, FDR = classC$FDR,
                            FDRsig = ifelse(classC$FDR < 0.05, TRUE, FALSE),
                            rank = NA, pairs = paste0(classC$signature, "_", classC$drug),
                            pset = "Meta Estimate", sig = NA, meta = TRUE))

df$pairs <- gsub("_", " and ", df$pairs)


###########################################################
# Plot Class C Biomarkers
###########################################################

png("DrugResponse/results/figures/robust_PSet/classC.png", width = 11, height = 12, res = 600, units = "in")
ggplot(df, aes(x = ci, y = pset)) + geom_linerange(aes(xmin = lower, xmax = upper)) + 
    geom_vline(xintercept = 0.5) + geom_vline(xintercept = c(0.4, 0.6), linetype = "dotted") + 
    geom_point(data = df, aes(shape = meta), size=8, fill = "#DB504A") +
    scale_shape_manual(values=c(23, 15), labels = c("Meta Estimate", "Concordance Index")) + 
    scale_y_discrete(limits=rev) +
    guides(shape = guide_legend(title = NULL)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0.02), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    facet_wrap(pairs ~ ., ncol = 3, nrow = 4) + theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          axis.text.x = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16), 
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          strip.text = element_text(size = 14),
          axis.title.y = element_text(size = 16)) +
    labs(x = "Concordance Index | Meta Estimate", y = "PSet", fill = "FDR Significance")
dev.off()