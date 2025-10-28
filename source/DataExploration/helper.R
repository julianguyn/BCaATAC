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

#' Compute BCa subtypes using genefu
#' 
#' @param rna dataframe. RNA-Seq counts matrix
#' @param meta dataframe. Gene metadata
#' @param model string. Subtyping classification model for genefu
#' Options: "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow"
#' @return dataframe of subtype scores per sample
score_bcasubtype <- function(rna, meta, model) {

    # format matrix
    rna <- t(rna) |> as.data.frame()

    # set do.mapping param
    do.mapping <- ifelse(model == "pam50", FALSE, TRUE)

    if (do.mapping == TRUE) {   # if not PAM50
        subset_meta <- meta[-which(is.na(meta$EntrezGene.ID)),]
        rna <- rna[,match(subset_meta$Ensembl, colnames(rna))]
        colnames(rna) <- subset_meta$EntrezGene.ID[match(subset_meta$Ensembl, colnames(rna))]

        # average across duplicates
        unique_genes <- unique(colnames(rna))
        rna <- as.data.frame(sapply(unique_genes, function(x) {
            rowMeans(rna[, names(rna) == x, drop = FALSE])
        }))
        meta <- subset_meta[!duplicated(subset_meta$EntrezGene.ID),]
        dimnames(meta)[[1]] <- meta$EntrezGene.ID
    } else {                    # if PAM50
        # match and get gene names
        colnames(rna) <- meta$Gene.Symbol[match(colnames(rna), meta$Ensembl)]
    }

    # get subtype predictions for PAM50
    SubtypePredictions <- molecular.subtyping(sbt.model = model, 
                                              data = rna,
                                              annot = meta, 
                                              do.mapping = do.mapping)

    # reutrn subtypes and subtyping probability
    res_df <- cbind(Subtype = as.character(SubtypePredictions$subtype), 
                    as.data.frame(SubtypePredictions$subtype.proba))
    return(res_df)
}

#' Get cell line subtype
#'
get_cell_subtype <- function() {

    anno <- read.csv("MetaData/Lupien/BCa_samples.csv")
    cells <- anno[anno$type == "cell.line",]
    cells <- cells[-which(cells$dup == 'T' & cells$dup_nerg == FALSE),]
    cells$filename <- gsub("_peaks.*", "", cells$filename)

    samples <- get_cells()
    samples$file <- gsub("\\.(?!$)", "-", samples$file, perl = TRUE)
    cells <- cells[cells$filename %in% samples$file,]
    cells$sample <- samples$sample[match(cells$filename, samples$file)]

    # set 600MPE as Luminal A: https://pmc.ncbi.nlm.nih.gov/articles/PMC5001206/
    cells$subtype[cells$sample == "600MPE"] <- "LumA"

    return(cells[,c(1:5)])
}
