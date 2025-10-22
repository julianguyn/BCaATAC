#' List of BCa-relevant drugs
#'
bca_drugs <- c("5-Fluorouracil",
               "Abemaciclib",
               "Alpelisib",
               "Anastrozole",
               "AZD-5305",
               "AZD-8205",
               "BMN-673",
               "Capecitabine",
               "Carboplatinum",
               "CDX-011",
               "CFI-400945",
               "CFI-402257",
               "Cisplatin",
               "CX-5461",
               "Cyclophosphamide",
               "Docetaxel",
               "Doxorubicin",
               "Elacestrant",
               "Epirubicin",
               "Eribulin",
               "Erlotinib",
               "Everolimus",
               "Exemestane",
               "Fluvastatin",
               "Fulvestrant",
               "Gemcitabine",
               "Inavolisib",
               "Ipatasertib",
               "Lapatinib",
               "Letrozole",
               "Neratinib",
               "Olaparib",
               "Paclitaxel",
               "Palbociclib",
               "Pembrolizumab",
               "Pertuzumab",
               "Ribociclib",
               "Selumetinib",
               "Tamoxifen",
               "Trastuzumab",
               "Tucatinib",
               "UNC0642")

#' Search for BCa drugs in PSets
#'
#' Searches for BCa drugs in pset drug sensitivity data.
#' @param sen dataframe. Drug sensitivity dataframe
#' @param label string. PSet label
#' @return A dataframe of meta analysis results.
#' 
drug_overlap <- function(sen, label) {
    tested_drugs <- rownames(sen)
    df <- data.frame(Drug = bca_drugs, Present = "Absent")

    for (i in 1:length(bca_drugs)) {
        drug = bca_drugs[i]
        if (drug %in% tested_drugs) {
            df$Present[i] <- "Present"
        }
    }
    overlap <- tested_drugs[tested_drugs %in% bca_drugs]
    print(paste("Number of overlapping drugs:", length(overlap)))
    df$PSet <- label
    return(df)
}

#' Creates dataframe for drug correlation and plotting
#'
#' @param pset1 dataframe. Drug sensitivity dataframe1
#' @param pset2 dataframe. Drug sensitivity dataframe2
#' @return dataframe object for correlation and plotting.
#' 
format_drug_pset <- function(pset1, pset2) {

    # get common samples
    common_cell <- intersect(colnames(pset1), colnames(pset2))
    common_drug <- intersect(rownames(pset1), rownames(pset2))
    pset1 <- pset1[match(common_drug, rownames(pset1)),match(common_cell, colnames(pset1))]
    pset2 <- pset2[match(common_drug, rownames(pset2)),match(common_cell, colnames(pset2))]

    # set labels and melt
    pset1$drug <- rownames(pset1)
    pset2$drug <- rownames(pset2)
    pset1 <- melt(pset1)
    pset2 <- melt(pset2)
    pset1$pairs <- paste0(pset1$variable, "_", pset1$drug)
    pset2$pairs <- paste0(pset2$variable, "_", pset2$drug)

    # # scatter plot of drug response difference for CTRP and GDSC 
    toPlot <- data.frame(
        pair = pset1$pairs, 
        pset1 = pset1$value, 
        pset2 = pset2$value[match(pset1$pairs, pset2$pairs)]
    )
    
    return(toPlot)
}

#' Load in RNA-Seq counts matrix from UHNBreast2 PSet
#' 
load_bca_RNA <- function() {

    # get RNA-seq matrix
    ubr2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
    ubr2 <- ubr2@molecularProfiles@ExperimentList$genes_counts
    rna <- ubr2@assays@data$expr
    colnames(rna) <- ubr2@colData$sampleid
    #rownames(rna) <- ubr2@rowRanges$gene_name

    # from map_sen()
    # missing: "HBL100" "HCC1008" 
    for (i in 1:length(colnames(rna))) {
        cell = colnames(rna)[i]
        if (cell %in% names(mapping_cells)) {colnames(rna)[i] <- unname(mapping_cells[cell])}
    }

    # keep only cell lines being used
    samples <- get_cells()
    rna <- rna[,colnames(rna) %in% samples$sample]
    df <- cbind(data.frame(Genes = rownames(rna), rna))
    write.table(df, file = "Signatures/data/bcacells_gene_counts.matrix",  quote = F, sep = "\t", col.names = T, row.names = F)

    return(rna)
}

#' Load in RNA-Seq counts matrix from other PSets
#' 
get_pset_rna <- function(filepath) {
    pset <- readRDS(filepath) |> updateObject()
    pset <- summarizeMolecularProfiles(pset, mDataType = "Kallisto_0.46.1.rnaseq.counts")
    rna <- pset@assays@data$expr
    #rownames(rna) <- rowData(pset)$gene_name

    # keep only cell lines being used
    samples <- get_cells()
    rna <- rna[,colnames(rna) %in% samples$sample]
    return(rna)
}

#' Correlate RNA-Seq expression of two psets
#' 
#' @param pset1 dataframe. Melted gene counts matrix of first pset
#' @param pset2 dataframe. Melted gene counts matrix of second pset
#' 
corr_pset_rna <- function(pset1, pset2) {

    # keep common cell lines
    ccls <- intersect(pset1$Var2, pset2$Var2)

    pset1 <- pset1[pset1$Var2 %in% ccls,]
    pset2 <- pset2[pset2$Var2 %in% ccls,]

    # check that order is the same
    pset1$pairs <- paste0(pset1$Var1, pset1$Var2)
    pset2$pairs <- paste0(pset2$Var1, pset2$Var2)
    table(pset1$pairs == pset2$pairs)

    # pearson correlation coefficient
    corr <- cor(pset1$value, pset2$value,  method = "spearman", use = "complete.obs")
    return(corr)
}
