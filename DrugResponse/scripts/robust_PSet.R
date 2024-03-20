setwd("C:/Users/julia/Documents/BCaATAC")

library(ggplot2)
library(ggh4x)
suppressMessages(library(meta))

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93")
pal2 <- c("#899DA4", "#BC4749")

# load in drug response associations
load("DrugResponse/results/data/indiv_PSet_CI.RData")
all_com <- rbind(ubr1_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com) #combine all drug response from all psets

sig_com$pairs <- paste0(sig_com$signature, "_", sig_com$drug)
all_com$pairs <- paste0(all_com$signature, "_", all_com$drug)


##########################
### CLASS B BIOMARKERS ###
##########################

# keep only signature-drug pairs that are significant in 2+ PSets
df <- all_com[which(all_com$pairs %in% unique(sig_com$pairs[duplicated(sig_com$pairs)])),]
df$sig <- ifelse(df$ci > 0.6 | df$ci < 0.4, "sig", "not sig")
df[df$FDRsig == FALSE,]$sig <- "not sig"
df$pset <- toupper(df$pset)


# specify sensitivity and resistance
df$type <- ifelse(df$ci > 0.5, "Sensitivity", "Resistance")
df$type <- factor(df$type, levels = c("Sensitivity", "Resistance"))


# plot all 
png("DrugResponse/results/figures/robust_PSet/two_psets_all.png", width = 16, height = 5, res = 600, units = "in")
ggplot(df, aes(x = pset, y = ci - 0.5, fill = ifelse(sig == "sig", type, "Not Significant"))) + geom_bar(stat="identity", color = "black") +
    facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0), labels = function(y) y + 0.5, oob = scales::squish) +
    scale_fill_manual(values = c(pal2, "#E8E1D9"), labels = c("Sensitivity", "Resistance", "Not Significant")) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "PSet") + 
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") + geom_hline(yintercept = c(0)) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# identify class B
df$classB <- FALSE

for (pair in unique(df$pairs)) {
    # for each pair, check if sig associations are all the same type (all sensitivity or all resistance)
    tmp <- df[df$pairs == pair,]
    ifelse(length(unique(tmp[tmp$sig == "sig",]$type)) == 1, df[df$pairs == pair,]$classB <- TRUE, df[df$pairs == pair,]$classB <- FALSE)
}

# remove non-significant associations and only plot classB

df <- df[df$classB == TRUE,]
df <- df[df$sig == "sig",]

# plot class B
png("DrugResponse/results/figures/robust_PSet/classB.png", width = 12, height = 4, res = 600, units = "in")
ggplot(df, aes(x = pset, y = ci - 0.5, fill = signature)) + geom_bar(stat="identity", color = "black") +
    facet_nested(~ factor(signature) + factor(drug), scales = "free_x") +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0), labels = function(y) y + 0.5, oob = scales::squish) +
    scale_fill_manual(values = pal) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "PSet") + 
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") + geom_hline(yintercept = c(0)) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


##########################
### CLASS C BIOMARKERS ###
##########################

# keep only signature-drug pairs that are present in at least 3 PSets
df <- all_com[which(all_com$pairs %in% names(table(all_com$pairs)[table(all_com$pairs) > 2])),]
df$sig <- ifelse(df$ci > 0.6 | df$ci < 0.4, "sig", "not sig")
df[df$FDRsig == FALSE,]$sig <- "not sig"
df$pset <- toupper(df$pset)

# data frame to hold meta estimates
estimates <- as.data.frame(matrix(data = NA, nrow = length(unique(df$pairs)), ncol = 9))
colnames(estimates) <- c("signature","drug", "TE", "seTE", "upper", "lower", "pval", "pval.Q", "I2")


for (i in 1:length(unique(df$pair))) {
  print(paste(i, "out of", length(unique(df$pair))))
  
  pair <- unique(df$pairs)[i]
  
  tmp <- df[which(df$pair == pair),]
  head(tmp)
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
  #pdf(fileName, width=10 , height=8, onefile=FALSE)
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
write.csv(estimates, file = "DrugResponse/results/tables/robust_PSet/meta_estimates.csv", quote = F, row.names = F)

# get Class C biomarkers
classC <- estimates[which(estimates$TE > 0.6 | estimates$TE < 0.4),]
classC <- classC[classC$FDR < 0.05,]

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

png("DrugResponse/results/figures/robust_PSet/classC.png", width = 10, height = 6, res = 600, units = "in")
ggplot(df, aes(x = ci, y = pset)) + geom_linerange(aes(xmin = lower, xmax = upper)) + 
    geom_vline(xintercept = 0.5) + geom_vline(xintercept = c(0.4, 0.6), linetype = "dotted") + 
    geom_point(data = df, aes(fill = FDRsig, shape = meta), size=5) +
    scale_shape_manual(values=c(18, 22)) + scale_y_discrete(limits=rev) +
    guides(shape = "none", fill = guide_legend(override.aes=list(shape=22, size = 8))) +
    scale_fill_manual(guide = guide_legend(reverse = FALSE), values = unname(pal), label = c("Significant", "Not Significant")) + 
    scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0.02), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    facet_wrap(pairs ~ .) + theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    labs(x = "Concordance Index (CI) | Meta Estimate", y = "PSet", fill = "FDR Significance")
dev.off()



# function to plot individual signatures
plot_signature <- function(signature) {
    # subset dataframe to keep only signature of interest
    toPlot <- df[df$signature == signature,]

    # new facet label names for facets
    labs <- unique(toPlot$pairs)
    supp.labs <- setNames(gsub(".*_", "", labs), labs)

    # create plot
    p <- ggplot(toPlot, aes(x = pset, y = ci - 0.5, fill = ifelse(sig == "sig", signature, "Not Significant"))) + geom_bar(stat="identity", color = "black") +
    facet_grid(~ factor(pairs), scales = "free_x", labeller = as_labeller(supp.labs)) +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0), labels = function(y) y + 0.5, oob = scales::squish) +
    scale_fill_manual(values = pal) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") + geom_hline(yintercept = c(0)) +
    theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))

    return(p)
}

# Signature 0
png("DrugResponse/results/figures/robust_PSet/Lucie/sig0.png", width = 16, height = 5, res = 600, units = "in")
plot_signature("Signature0")
dev.off()

# Signature 1
png("DrugResponse/results/figures/robust_PSet/Lucie/sig1.png", width = 16, height = 5, res = 600, units = "in")
plot_signature("signature1")
dev.off()

# Signature 2
png("DrugResponse/results/figures/robust_PSet/Lucie/sig2.png", width = 16, height = 5, res = 600, units = "in")
plot_signature("signature2")
dev.off()

# Signature 3
png("DrugResponse/results/figures/robust_PSet/Lucie/sig3.png", width = 16, height = 5, res = 600, units = "in")
plot_signature("signature3")
dev.off()

# Signature 4
png("DrugResponse/results/figures/robust_PSet/Lucie/sig4.png", width = 16, height = 5, res = 600, units = "in")
plot_signature("signature4")
dev.off()

# Signature 5
png("DrugResponse/results/figures/robust_PSet/Lucie/sig5.png", width = 16, height = 5, res = 600, units = "in")
plot_signature("signature5")
dev.off()
