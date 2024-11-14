setwd("C:/Users/julia/Documents/BCaATAC")

library(readxl)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# read in drug response data
response <- as.data.frame(read_excel("PacResistantPDX/data/PDX_Data_Xeva.xlsx", sheet = 1))

# keep only mRECIST
response <- response[,colnames(response) %in% c("patient_id", "mRECIST", "drug")]


# read in signature scores
scores <- as.data.frame(t(read.table("PacResistantPDX/data/zscore.tsv")))
rownames(scores) <- toupper(gsub("_", "", gsub("X", "", rownames(scores))))

# subset for just paclitaxel and for samples with signature scores
pac <- response[response$drug == "TAXOL",] # removed TAXOL CONTROL?
pac <- pac[pac$patient_id %in% rownames(scores),]
pac <- pac[-which(pac$mRECIST == "NA"),]
pac$mRECIST <- factor(pac$mRECIST, levels = c("CR", "PR", "SD", "PD"))

# loop through each signature and create the waterfall plot
for (i in 1:ncol(scores)) {
    
    # save signature name
    signature <- colnames(scores)[i]

    # reset signature score
    pac$signature_score <- NA

    # get signature score
    for (j in 1:nrow(pac)) { pac$signature_score[j] <- scores[rownames(scores) == pac$patient_id[j],i] }

    # formating for plotting
    pac <- pac[order(pac$signature_score, decreasing = T),]
    pac$rank <- 1:nrow(pac)


    # get min max values for y axis
    val <- round(max(abs(min(pac$signature_score)), max(pac$signature_score))) + 5

    # legend labels and colors
    labs <- c('CR' = 'Complete\nResponse', 'PD' = 'Progressive\nDisease', 'SD' = 'Stable\nDisease','PR' = 'Partial\nResponse')
    pal <- c("#136F63", "#9DCBBA", "#FFADA1", "#B02E0C")

    # waterfall plot
    png(paste0("PacResistantPDX/results/figures/", signature ,".png"), width=175, height=150, units='mm', res = 600, pointsize=80)
    print(ggplot(pac, aes(x = rank, y = signature_score, fill = mRECIST)) + 
        geom_col(color = "black") + geom_hline(yintercept = 0) +
        scale_x_discrete(limits = 1:length(pac$rank), labels = pac$patient_id) +
        scale_fill_manual(values = pal, labels = labs) +
        theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #, axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylim(c(-val, val)) + labs(x = "PDX Model", y = "Signature Score", fill = "mRECIST", title = signature) )
    dev.off()
}



####### Try with best average response

setwd("C:/Users/julia/Documents/BCaATAC")

library(readxl)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# read in drug response data
response <- as.data.frame(read_excel("PacResistantPDX/data/PDX_Data_Xeva.xlsx", sheet = 1))

# keep only best.average.response
response <- response[,colnames(response) %in% c("patient_id", "best.average.response", "drug")]


# read in signature scores
scores <- as.data.frame(t(read.table("PacResistantPDX/data/zscore.tsv")))
rownames(scores) <- toupper(gsub("_", "", gsub("X", "", rownames(scores))))

# subset for just paclitaxel and for samples with signature scores
pac <- response[response$drug == "TAXOL",] # removed TAXOL CONTROL?
pac <- pac[pac$patient_id %in% rownames(scores),]
pac <- pac[-which(pac$best.average.response == "NA"),]
pac$best.average.response <- as.numeric(pac$best.average.response)

# loop through each signature and create the waterfall plot
for (i in 1:ncol(scores)) {
    
    # save signature name
    signature <- colnames(scores)[i]

    # reset signature score
    pac$signature_score <- NA

    # get signature score
    for (j in 1:nrow(pac)) { pac$signature_score[j] <- scores[rownames(scores) == pac$patient_id[j],i] }
    pac$sig <- NA
    pac$sig <- ifelse(pac$signature_score > 0, "Positive", "Negative")
    pac$sig <- factor(pac$sig, levels = c("Positive", "Negative"))

    # formating for plotting
    pac <- pac[order(pac$best.average.response, decreasing = T),]
    pac$rank <- 1:nrow(pac)

    # get min max values for y axis
    val <- round(max(abs(min(pac$best.average.response)), max(pac$best.average.response))) + 5

    # waterfall plot 1
    png(paste0("PacResistantPDX/results/figures/BAR-1/", signature ,".png"), width=175, height=125, units='mm', res = 600, pointsize=80)
    print(ggplot(pac, aes(x = rank, y = best.average.response, fill = signature_score)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_gradientn(colours = brewer.pal(9, "Blues")) +
        theme_classic() + theme(legend.key.size = unit(0.6, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylim(c(-val, val)) + labs(x = "PDX Model", y = "Best Average Response", fill = "Signature\nScore", title = signature) )
    dev.off()


    # waterfall plot 2
    png(paste0("PacResistantPDX/results/figures/BAR-2/", signature ,".png"), width=175, height=125, units='mm', res = 600, pointsize=80)
    print(ggplot(pac, aes(x = rank, y = best.average.response, fill = sig)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_manual(values = c("#136F63", "#FFADA1")) +
        theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylim(c(-val, val)) + labs(x = "PDX Model", y = "Best Average Response", fill = "Signature\nScore", title = signature) )
    dev.off()
}

