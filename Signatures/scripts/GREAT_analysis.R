# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(rGREAT)
    library(wesanderson)
    library(ggplot2)
})

set.seed(123)

###########################################################
# Load in data
###########################################################

# read in background
bg <- as.data.frame(fread("Signatures/results/data/beds/background.bed"))

for (i in 1:6) {

    print("Starting Signature")
    # read in bed file
    gr <- as.data.frame(fread(paste0("Signature", i, ".bed")))
    #gr <- as.data.frame(fread(paste0("Signatures/results/data/beds/Signature", i, ".bed")))

    # submit job
    print("Submitting Job")
    job = submitGreatJob(gr, bg, species = "hg38")
    print("Job complete")
    tbl = getEnrichmentTables(job)

    # save results
    mf <- as.data.frame(tbl[1])[,c(1:2, 11:13)]
    mf$lab <- "Molecular Feature"
    bp <- as.data.frame(tbl[2])[,c(1:2, 11:13)]
    bp$lab <- "Biological Process"
    cc <- as.data.frame(tbl[3])[,c(1:2, 11:13)]
    cc$lab <- "Cellular Component"

    # add column names
    cols <- c("GO.ID", "GO.Name", "Total_Genes_Annotated", "Raw_PValue", "Adjp_BH", "Label")

    colnames(mf) <- cols
    colnames(bp) <- cols
    colnames(cc) <- cols

    # combine all results
    res <- rbind(bp, mf, cc)

    # filter results
    res <- res[res$Adjp_BH < 0.05,]
    res <- res[order(res$Adjp_BH),]

    print("Saving file")
    write.table(res, file = paste0("Signature", i, "_GREAT.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
}


plot_GREAT <- function(i) {
    df <- as.data.frame(fread(paste0("Signatures/results/data/GREAT/Signature", i, "_GREAT.tsv")))
    df$Total_Genes_Annotated <- as.numeric(df$Total_Genes_Annotated)
    df$Adjp_BH <- as.numeric(df$Adjp_BH)

    # keep only top 30
    toPlot <- df[1:30,]
    toPlot <- toPlot[order(toPlot$Adjp_BH, decreasing = TRUE),]
    toPlot$GO.Name <- factor(toPlot$GO.Name, levels=unique(toPlot$GO.Name))

    p <- ggplot(toPlot, aes(x = GO.Name, y = -log(Adjp_BH), size = Total_Genes_Annotated, color = Label), shape = 19) +
        geom_point() + #scale_size(range = c(0.1, 3)) +
        coord_cartesian(clip = "off") + coord_flip()  +
        guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=19, size = 4))) +
        #guides(shape = guide_legend(override.aes = list(size = 1)), color = guide_legend(override.aes = list(shape = 19, size = 1))) +
        scale_color_manual(values = c("#574B60", "#CBC9AD", "#5B96AF"), labels = c("BP", "CC", "MF")) +
        theme_classic() + 
        theme(legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            legend.position = c(0.76, 0.2), #plot.margin=unit(c(0.6,0,0,1),"cm"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + 
        labs(x = "GO Term", y = "-log(FDR)", size = "N", color = "Ontology")
    return(p)
}
png(paste0("Signatures/results/figures/Sig1_GREAT.png"), width = 7, height = 7, res = 600, units = "in")
plot_GREAT(1)
dev.off()

png(paste0("Signatures/results/figures/Sig2_GREAT.png"), width = 7, height = 7, res = 600, units = "in")
plot_GREAT(2)
dev.off()

png(paste0("Signatures/results/figures/Sig3_GREAT.png"), width = 7, height = 7, res = 600, units = "in")
plot_GREAT(3)
dev.off()

png(paste0("Signatures/results/figures/Sig4_GREAT.png"), width = 7, height = 7, res = 600, units = "in")
plot_GREAT(4)
dev.off()

png(paste0("Signatures/results/figures/Sig5_GREAT.png"), width = 7, height = 7, res = 600, units = "in")
plot_GREAT(5)
dev.off()

png(paste0("Signatures/results/figures/Sig6_GREAT.png"), width = 7, height = 7, res = 600, units = "in")
plot_GREAT(6)
dev.off()


