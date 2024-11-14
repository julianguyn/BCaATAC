setwd("C:/Users/julia/Documents/BCaATAC")

library(ggplot2)
library(RColorBrewer)
library(dplyr)

# read in drug response data
load("DrugResponse/results/data/sensitivity_data.RData")
colnames(ubr2_sen) <- gsub("-", "", colnames(ubr2_sen))

# read in signature scores
scores <- as.data.frame(t(read.table("PacResistantPDX/data/cells_zscore.tsv")))
rownames(scores)[rownames(scores) == "MDA_MB_175_VII_2_S5"] = "MDAMB175VII"
rownames(scores)[rownames(scores) == "LY2_1_S2"] = "MCF7/LY2"
rownames(scores)[rownames(scores) == "MDA361_2_S16"] = "MDAMB361"
rownames(scores)[rownames(scores) == "MB157."] = "MDAMB157"
rownames(scores)[rownames(scores) == "HS578."] = "Hs 578T"

rownames(scores) <- gsub("\\.", "", gsub("_.*", "", rownames(scores)))

##################
### UHNBREAST2 ###
##################

ubr2 <- scores[rownames(scores) %in% colnames(ubr2_sen),] #41

# format drug response 
ubr2_sen <- as.data.frame(t(ubr2_sen))
ubr2_sen <- ubr2_sen[rownames(ubr2_sen) %in% rownames(scores),]


# loop through each signature and create the waterfall plot
for (i in 1:ncol(scores)) {
    
    # save signature name
    signature <- colnames(scores)[i]

    # reset signature score
    ubr2_sen$signature_score <- NA

    # get signature score
    for (j in 1:nrow(ubr2_sen)) { ubr2_sen$signature_score[j] <- scores[rownames(scores) == rownames(ubr2_sen)[j],i] }

    # formating for plotting
    ubr2_sen <- ubr2_sen[order(ubr2_sen$Paclitaxel, decreasing = T),]
    ubr2_sen$rank <- 1:nrow(ubr2_sen)

    # waterfall plot
    png(paste0("PacResistantPDX/results/figures/cells/", signature ,".png"), width=175, height=125, units='mm', res = 600, pointsize=80)
    print(ggplot(ubr2_sen, aes(x = rank, y = (Paclitaxel/100) - 0.5, fill = signature_score)) + 
        scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
        geom_bar(stat = "identity", color = "black") + scale_y_continuous(limits = c(-0.5, 0.5), labels = function(y) y + 0.5) + geom_hline(yintercept = 0) + 
        theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(x = "Cell Line", y = "AAC", fill = "Signature Score", title = signature) )
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

