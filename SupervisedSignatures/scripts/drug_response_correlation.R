# Script to identify drugs of interest for supervised signature extraction

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(plyr)
    library(ComplexHeatmap)
    library(circlize)
    library(reshape2)
    library(ggplot2)
    library(ggh4x)
    library(ggpubr)
    library(grid)
    library(gridExtra)
})


### ===== Get drugs of interest ===== ###

# read in PDX drug response data
pdx_response <- as.data.frame(read_excel("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", sheet = 1))
pdx_response[pdx_response == "NA"] <- NA
pdx_response[pdx_response == "TAXOL"] <- "PACLITAXEL"
pdx_response[pdx_response == "CARBOPLATIN"] <- "CARBOPLATINUM"
pdx_response <- na.omit(pdx_response)

# remove CONTROLs and H2O samples (left with n=23 drugs)
drugs_of_interest <- unique(pdx_response$drug[-which(pdx_response$drug %in% c("H2O", pdx_response$drug[grep("CONTROL", pdx_response$drug)]))])
pdx_response <- pdx_response[pdx_response$drug %in% drugs_of_interest,]

# keep only drugs with at least 10 samples (left with n=11 drugs)
drugs_of_interest <- names(table(pdx_response$drug)[table(pdx_response$drug) > 10])
pdx_response <- pdx_response[pdx_response$drug %in% drugs_of_interest,]


### ===== Get cell line drug response for drugs of interest ===== ### 

# load in cell line drug response data
load("DrugResponse/results/data/sensitivity_data.RData")

# load in and extract needed cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# keep only drugs of interest
ubr1_sen <- ubr1_sen[toupper(rownames(ubr1_sen)) %in% drugs_of_interest,colnames(ubr1_sen) %in% samples$sample] #taxol only, 43 cells
ubr2_sen <- ubr2_sen[toupper(rownames(ubr2_sen)) %in% drugs_of_interest,colnames(ubr2_sen) %in% samples$sample] #taxol, carboplatin, & eribulin, 42 cells
gray_sen <- gray_sen[toupper(rownames(gray_sen)) %in% drugs_of_interest,colnames(gray_sen) %in% samples$sample] #taxol & carboplatin, 41 cells
gcsi_sen <- gcsi_sen[toupper(rownames(gcsi_sen)) %in% drugs_of_interest,colnames(gcsi_sen) %in% samples$sample] #taxol only, 25 cells
gdsc_sen <- gdsc_sen[toupper(rownames(gdsc_sen)) %in% drugs_of_interest,colnames(gdsc_sen) %in% samples$sample] #taxol only, 34 cells
ctrp_sen <- ctrp_sen[toupper(rownames(ctrp_sen)) %in% drugs_of_interest,colnames(ctrp_sen) %in% samples$sample] #taxol & carboplatin, 32 cells
ccle_sen <- ccle_sen[toupper(rownames(ccle_sen)) %in% drugs_of_interest,colnames(ccle_sen) %in% samples$sample] #taxol only, 23 cells


### ===== Correlate drug response for common cell lines ===== ### 

# set palette for plotting
col_fun <- colorRamp2(c(-1, 0, 1), c("#A85751", "white", "#66999B"))

# function to create heatmap
plot_heatmap <- function(mat, legend_label) {
    plot <- Heatmap(mat,
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        heatmap_legend_param = list(
            title = legend_label,
            color_bar = "continuous"
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(round(mat[i, j], 2), x, y, gp = gpar(fontsize = 8))
        }
    )
    return(plot)
}

# get all paclitaxel response
pac <- rbind.fill(
    ubr1_sen[rownames(ubr1_sen) == "Paclitaxel",],
    ubr2_sen[rownames(ubr2_sen) == "Paclitaxel",],
    gray_sen[rownames(gray_sen) == "Paclitaxel",],
    gcsi_sen[rownames(gcsi_sen) == "Paclitaxel",],
    gdsc_sen[rownames(gdsc_sen) == "Paclitaxel",],
    ccle_sen[rownames(ccle_sen) == "Paclitaxel",],
    ctrp_sen[rownames(ctrp_sen) == "Paclitaxel",]
) |> t() |> as.data.frame()
colnames(pac) <- c("UBR1", "UBR2", "GRAY", "gCSI", "GDSC", "CCLE", "CTRP")

# create correlation matrix
corr_mat <- cor(pac, use = "pairwise.complete.obs", method = "pearson")

# correlation heatmap for paclitaxel
png("SupervisedSignatures/results/figures/dr_corr/pcc-paclitaxel.png", width=125, height=100, units='mm', res = 600, pointsize=80)
plot_heatmap(corr_mat, "Pearson's")
dev.off()


# get all carboplatin response
car <- rbind.fill(
    ubr2_sen[rownames(ubr2_sen) == "Carboplatinum",],
    gray_sen[rownames(gray_sen) == "Carboplatinum",],
    ctrp_sen[rownames(ctrp_sen) == "Carboplatinum",]
) |> t()
colnames(car) <- c("UBR2", "GRAY", "CTRP")

# create correlation matrix
corr_mat <- cor(car, use = "pairwise.complete.obs", method = "pearson")

# correlation heatmap for carboplatinum
png("SupervisedSignatures/results/figures/dr_corr/pcc-carboplatinum.png", width=85, height=65, units='mm', res = 600, pointsize=80)
plot_heatmap(corr_mat, "Pearson's")
dev.off()


### ===== Concordance of binarized drug response values in cell lines ===== ### 

# function to binarize drug repsonse
binarize_dr <- function(sensitivity_data) {
    return(as.data.frame(ifelse(sensitivity_data >= 0.5, 1, 0)))
}

# binarize drug repsonse matrices
pac_bin <- binarize_dr(pac)
car_bin <- binarize_dr(car)


# function to compute jaccard similarity 
jaccard <- function(a, b) {
    # remove cell lines with NAs
    keep <- !is.na(a) & !is.na(b)
    a <- a[keep]
    b <- b[keep]

    intersection = sum(a == b)
    union = length(a)
    return (intersection/union)
}

# function to compute a matrix of jaccard similarity
jaccard_mat <- function(drug_response_mat) {

    # initialize matrix
    mat <- matrix(NA, nrow = ncol(drug_response_mat), ncol = ncol(drug_response_mat), 
                dimnames = list(names(drug_response_mat), names(drug_response_mat)))

    # compute pairwise jaccard similarity
    combn(seq_len(ncol(drug_response_mat)), 2, 
        function(idx) {
            i <- idx[1]
            j <- idx[2]
            
            # compute jaccard and fill in symmetric matrix
            mat[i, j] <<- mat[j, i] <<- jaccard(drug_response_mat[[i]], drug_response_mat[[j]])
        }
    )

    # set diagonal as 1
    diag(mat) <- 1

    return(mat)
}


# create jaccarad matrices
jacc_pac <- jaccard_mat(pac_bin)
jacc_car <- jaccard_mat(car_bin)

# correlation heatmap for paclitaxel
png("SupervisedSignatures/results/figures/dr_corr/jacc-paclitaxel.png", width=125, height=100, units='mm', res = 600, pointsize=80)
plot_heatmap(jacc_pac,"Jaccard\nSimilarity")
dev.off()

# correlation heatmap for paclitaxel
png("SupervisedSignatures/results/figures/dr_corr/jacc-carboplatinum.png", width=85, height=65, units='mm', res = 600, pointsize=80)
plot_heatmap(jacc_car,"Jaccard\nSimilarity")
dev.off()


### ===== Compute overall paclitaxel responses ===== ### 

# compute average drug response (AAC)
pac$response <- rowMeans(pac, na.rm = TRUE)

# compute average drug response (binarized response)
pac_bin$avg <- rowMeans(pac_bin, na.rm = TRUE)
pac_bin$response <- ifelse(pac_bin$avg > 0.5, "Sensitive", ifelse(pac_bin$avg < 0.5, "Not Sensitive", NA))


# save binarized paclitaxel drug response
save(pac, file = "SupervisedSignatures/results/data/pac-response.RData")
save(pac_bin, file = "SupervisedSignatures/results/data/pac-binarized.RData")


### ===== Plot heatmap of paclitaxel response (AAC) ===== ### 

# get cell line clustering order (from signature scoring)
order <- c("HCC70", "MDA-MB-157", "CAL-85-1", "HCC38", "MDA-MB-468", "HBL-100", "HCC1187", 
           "CAL-51", "DU4475", "MCF10A", "SK-BR-7", "HDQ-P1", "JIMT-1", "SUM149PT", 
           "HCC1806", "HCC1143", "BT-20", "HCC1937", "CAL-120", "HCC1395", "MDA-MB-436", 
           "BT-549", "Hs 578T", "MDA-MB-231", "SUM159PT", "MDA-MB-175-VII", "MDA-MB-361", 
           "MX-1", "HCC3153", "MCF-7", "KPL-1", "MCF-7/LY2", "HCC1954", "HCC1008", "SK-BR-3", 
           "CAL-148", "AU565", "BT-474", "CAMA-1", "EFM-192A", "HCC202", "UACC-812", 
           "SK-BR-5", "600MPE", "EFM-19", "ZR-75-1", "T-47D", "HCC1428", "SUM52PE")
order <- order[order %in% rownames(pac_bin)]

# set colours for plotting
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")

# read in subtype annotation
meta <- fread("DataExploration/data/bcacells_annotation.tsv")
anno_plot <- data.frame(Cell = meta$sample, Subtype = meta$subtype)
anno_plot <- anno_plot[match(order, anno_plot$Cell),]
anno_plot$Label <- "Subtype"
anno_plot$Cell <- factor(anno_plot$Cell, levels = order)

# format dataframe for plotting
pac$Cell <- rownames(pac)
toPlot <- melt(pac[,-which(colnames(pac) == "response")])
toPlot$Cell <- factor(toPlot$Cell, levels = order)

# plot individual response per PSet per cell line
p1 <- ggplot(toPlot, aes(x = variable, y = Cell, fill = value)) + 
  geom_tile(color = "white") + theme_void() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradientn(colors = c("#BC4749", "#A69AD8", "#7294D4"), na.value = "#D1D7DD") +        
  theme(
        axis.text.x = element_text(vjust = 0.5),
        axis.title.x = element_text(),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    ) +
  labs(x = "\nPSet", y = "Cell Line", fill = "Response\n(AAC)")

# plot overall paclitaxel response
toPlot <- data.frame(Cell = rownames(pac), Response = pac$response)
toPlot$Label <- "Overall"
toPlot$Cell <- factor(toPlot$Cell, levels = order)

p2 <- ggplot(toPlot, aes(x = Label, y = Cell, fill = Response)) + 
  geom_tile(color = "white") + theme_void() +
  geom_text(aes(label = round(Response, 2)), color = "black", size = 3) +
  scale_fill_gradientn(colors = c("#BC4749", "#A69AD8", "#7294D4"), na.value = "#D1D7DD") + 
  theme(
        axis.text.x = element_text(vjust = 0.5),
        axis.title.x = element_text(),
        strip.text.x = element_text(),
        legend.position = "none"
    ) +
  labs(x = "\n", fill = "")

# plot subtypes
p3 <- ggplot(anno_plot, aes(x = Label, y = Cell, fill = Subtype)) + 
  geom_tile(color = "white") + theme_void() +
  scale_fill_manual(values = subtype_pal) +
  theme(
        axis.text.x = element_text(vjust = 0.5),
        axis.title.x = element_text(),
        strip.text.x = element_text()
    ) +
  labs(x = "\n", fill = "Subtype")

# extract legends
l1 <- as_ggplot(get_legend(p1))
l3 <- as_ggplot(get_legend(p3))
p1 <- p1+theme(legend.position = "none")
p3 <- p3+theme(legend.position = "none")

png("SupervisedSignatures/results/figures/pac_sen_AAC_heatmap.png", width=200, height=200, units='mm', res = 600, pointsize=80)
grid.arrange(p1, p2, p3, l1, l3, ncol = 30, nrow = 3,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,NA,NA,NA,NA,NA),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,4,4,4,4,4),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,5,5,5,5,5),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,NA,NA,NA,NA,NA)))
dev.off()



### ===== Plot heatmap of binarized paclitaxel response ===== ### 

# format dataframe for plotting
pac_bin$Cell <- rownames(pac_bin)
toPlot <- melt(pac_bin[,-which(colnames(pac_bin) %in% c("avg", "response"))])
toPlot$value <- factor(toPlot$value, levels = c("1", "0", "NA"))
toPlot$Cell <- factor(toPlot$Cell, levels = order)

# plot individual response per PSet per cell line
p1 <- ggplot(toPlot, aes(x = variable, y = Cell, fill = value)) + 
  geom_tile(color = "white") + theme_void() +
  scale_fill_manual(values = c("#7294D4", "#BC4749"), labels = c("Sensitive", "Not Sensitive"), na.value = "#D1D7DD") +
  theme(
        axis.text.x = element_text(vjust = 0.5),
        axis.title.x = element_text(),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    ) +
  labs(x = "\nPSet", y = "Cell Line", fill = "")

# plot overall paclitaxel response
toPlot <- data.frame(Cell = rownames(pac_bin), Response = pac_bin$response)
toPlot$Response <- factor(toPlot$Response, levels = c("Sensitive", "Not Sensitive"))
toPlot$Label <- "Overall"
toPlot$Cell <- factor(toPlot$Cell, levels = order)

p2 <- ggplot(toPlot, aes(x = Label, y = Cell, fill = Response)) + 
  geom_tile(color = "white") + theme_void() +
  scale_fill_manual(values = c("#7294D4", "#BC4749"), na.value = "#D1D7DD") +
  theme(
        axis.text.x = element_text(vjust = 0.5),
        axis.title.x = element_text(),
        strip.text.x = element_text(),
        legend.position = "none"
    ) +
  labs(x = "\n", fill = "")

# extract legends
l1 <- as_ggplot(get_legend(p1))
p1 <- p1+theme(legend.position = "none")

png("SupervisedSignatures/results/figures/pac_sen_bin_heatmap.png", width=200, height=200, units='mm', res = 600, pointsize=80)
grid.arrange(p1, p2, p3, l1, l3, ncol = 30, nrow = 3,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,NA,NA,NA,NA,NA),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,4,4,4,4,4),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,5,5,5,5,5),
                          c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,NA,NA,NA,NA,NA)))
dev.off()