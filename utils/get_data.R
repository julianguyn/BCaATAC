#' Get unique 49 BCa cell lines for analysis
#'
#' For duplicated cell lines, keep only those from Nergiz
#' @return Subsetted dataframe of 49 sample filenames.
#' 
get_cells <- function() {

    # load in BCa cell lines
    samples <- read.csv("metadata/cl_samples.csv")

    # get duplicates
    dup <- samples$sample[duplicated(samples$sample)]
    tmp <- samples[samples$sample %in% dup,]

    # keep Nergiz dups
    samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),])

    return(samples)
}

#' Get ARCHE scores for samples
#'
#' @param sample string. "cells" or "pdxs"
#' @param arche string. "k20", "k50", or "all"
#' @return ARCHE scores for requested samples.
#' 
get_arche_scores <- function(sample, arche, meta) {

    # get filepath
    filepath <- switch(
        sample,
        cells = "data/rawdata/ccls/cells_",
        pdxs = "data/rawdata/pdx/PDXs_"
    )

    # get filename
    filename <- switch(
        arche,
        k20 = "20k.Zscore.txt",
        k50 = "50k.Zscore.txt",
        all = "all.Zscore.txt"
    )

    file <- paste0(filepath, filename)

    # load in arche scores
    scores <- read.table(file)
    rownames(scores) <- paste0("ARCHE", 1:6)

    # standardize sample names of cell lines
    colnames(scores) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(scores), perl = TRUE))
    scores <- scores[, colnames(scores) %in% meta$filename]
    colnames(scores) <- meta$sampleid[match(colnames(scores), meta$filename)]
    scores <- scores[, order(colnames(scores))]
    if ("104987" %in% colnames(scores)) scores <- scores[, colnames(scores) != "104987"]
    if (sample == "pdxs") colnames(scores) <- map_pdx(colnames(scores))

    return(scores)
}

#' Extract AAC for drug of interest
#'
#' Helper function for plot_CCLs().
#' @param pset dataframe. PSet drug sensitivity from get_drugsen()
#' @param drugs string. Drug of interest.
#' @param df dataframe. Dataframe to bind results to.
#' @param label string. PSet name.
#' @return df with drug AAC binded (if present in PSet).
#' 
get_doiAAC <- function(pset, drug, df, label) {
    if (drug %in% rownames(pset)) {
        # save drug response
        res <- reshape2::melt(pset[rownames(pset) == drug,]) |> suppressMessages()
        colnames(res) <- c("Sample", "AAC")
        res$PSet <- label
        # merge results with compiled dataframe
        df <- rbind(df, res)
        return(df)
    } else {return(df)}
}

#' Get CCL gene counts matrices
#' 
#' @param pset string. Name of PSet to load
#' @param gene.symbol boolean. True is rownames should be Gene Symbols
#' 
get_pset_rna <- function(pset, gene.symbol = FALSE) {

    pset_options <- c("UBR2", "GRAY", "gCSI", "CCLE")
    if (!pset %in% pset_options) {
        message("Option not valid. Pick from: UBR2, GRAY, gCSI, CCLE")
        return(NA)
    }

    filename <- paste0("data/procdata/CCLs/rna/", pset, "_RNA.tsv")
    df <- fread(filename, data.table = FALSE)
    df$Genes <- gsub("\\..*", "", df$Genes)

    if (gene.symbol == TRUE) {
        c_meta <- read.table("data/procdata/CCLs/rna/UBR2_RNA_meta.tsv", header = TRUE)
        df$Genes <- c_meta$Gene.Symbol[match(df$Genes, c_meta$Ensembl)]
    }

    df <- df %>%
        dplyr::group_by(Genes) %>%
        dplyr::summarise(across(everything(), mean)) %>%
        as.data.frame()
    rownames(df) <- df$Genes
    df$Genes <- NULL
    return(df)
}

#' Get cell line subtype
#'
get_cell_subtype <- function() {

    anno <- read.csv("metadata/BCa_samples.csv") # from lupien lab
    cells <- anno[anno$type == "cell.line",]
    cells <- cells[-which(cells$dup == 'T' & cells$dup_nerg == FALSE),]
    cells$filename <- gsub("_peaks.*", "", cells$filename)

    samples <- get_cells()
    samples$file <- gsub("\\.(?!$)", "-", samples$file, perl = TRUE)
    cells <- cells[cells$filename %in% samples$file,]
    cells$sample <- samples$sample[match(cells$filename, samples$file)]

    # set 600MPE as Luminal A: https://pmc.ncbi.nlm.nih.gov/articles/PMC5001206/
    # https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3
    cells$subtype[cells$sample == "600MPE"] <- "LumA"

    return(cells[,c(1:5)])
}

#' Get TDXd cell line IC50 response data
#' 
get_tdxd <- function() {
    tdxd <- read_excel("data/rawdata/psets/TDXd Cell Line Response Data.xlsx", sheet = 1) |>
        as.data.frame()
    tdxd <- tdxd[,colnames(tdxd) %in% c("Cell Line", "Average IC50...4", "Average IC50...9")]
    colnames(tdxd) <- c("Cell.Line", "Avg.IC50", "Avg.IC50.Treps")

    # remove no response values
    tdxd$Avg.IC50[tdxd$Avg.IC50 == "#DIV/0!"] <- NA
    tdxd$Avg.IC50.Treps[tdxd$Avg.IC50.Treps == "#DIV/0!"] <- NA
    return(tdxd)
}

#' Get ARCHE assignments from NMF output
#' 
#' Helper function to speed up loading in ARCHE matrix
#' @param mat boolean. If TRUE, return just the matrix itself
#' 
get_arche_tcga <- function(mat = FALSE) {

    # read in meta data
    meta <- read.csv("metadata/procdata/TCGA_subtype_label.csv")

    # load in matrix file from NMF

    atac <- read.table("data/procdata/NMF/ATAC_heatmap_rank6.png.order.matrix", header = T)
    
    if (mat == TRUE) {
        # if mat == TRUE, return just the matrix
        t_meta <- read.csv("data/rawdata/tcga/TCGA_sourcefiles.csv")
        colnames(atac) <- t_meta$Sample.Name[match(gsub("X", "", colnames(atac)), t_meta$ATAC.Seq.File.Name)]
        rownames(atac) <- paste0("ARCHE", 1:6)
        return(atac)

    } else {

        atac$ARCHE <- paste0("ARCHE", 1:6)
        atac$ARCHE <- factor(atac$ARCHE, levels = paste0("ARCHE",6:1))
        mat <- reshape2::melt(atac)
        mat$variable <- gsub("X", "", mat$variable)

        # save signature assignment
        mat$signature_assign <- ""
        for (sample in mat$variable) {
            tmp <- mat[mat$variable == sample,]
            mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$ARCHE)
        }
        #mat <- mat[,colnames(mat) %in% c("variable", "signature_assign")] |> unique()

        # save subtype
        mat$subtype <- meta$bca_subtype[match(gsub("X", "", mat$variable), gsub("-", "\\.", meta$File.Name))]

        mat$rank <- 1:nrow(mat)
        return(mat)
    }
}

#' Get TCGA RNA-Seq gene counts matrix
#' 
get_tcga_rna <- function() {
    counts <- t(fread("data/procdata/TCGA/TCGA_BRCA_gene_counts.matrix")) |> as.data.frame()
    colnames(counts) <- counts[1,]
    counts <- counts[-c(1),]
    genes <- rownames(counts)
    counts <- sapply(counts, as.double)
    rownames(counts) <- genes
    return(counts)
}

#' Get PDX ARCHE scores
#'
#' @return ARCHE scores fataframe with standardized PDX IDs.
#' 
get_arche_pdx <- function() {
    scores <- as.data.frame(t(read.table("data/rawdata/pdx/bca_sign.Zscore.txt")))
    colnames(scores) <- paste0("ARCHE", 1:6)
    rownames(scores) <- map_pdx(rownames(scores))
    return(scores)
}

#' Load in XEVA drug response spreadsheets
#' 
get_xeva <- function(type) {

    file <- switch(
        type,
        limited = "data/rawdata/pdx/PDX_Response_Sept2025.csv",
        full = "data/rawdata/pdx/PDX_BR_BAR_grouped.csv"
    )

    # read in file
    xeva <- read.csv(file)
    

    # create patient.id and standardize
    xeva$patient.id <- map_pdx(xeva$PDX_ID)

    # create instances for NOTCH01 passages
    notch4 <- notch6 <- xeva[xeva$patient.id %in% c("NOTCH01"),]
    notch4$patient.id <- "NOTCH01P4"
    notch6$patient.id <- "NOTCH01P6"
    xeva <- xeva[-which(xeva$patient.id %in% c("NOTCH01")),]
    xeva <- rbind(xeva, notch4, notch6)

    # create model_group
    xeva$model_group <- paste(xeva$patient.id, xeva$drug, sep = "_")

    return(xeva)
}


#' Get CpG sites from 450k array in each ARCHE
#' 
#' @param ann450k anno 
#' ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' 

# helper function to find cpgs overlapping in arches
anno_cpgs <- function(arche) {

    # gr for arche
    bed <- fread(paste0("data/procdata/ARCHEs/beds/", arche, "_20k.bed"))
    gr <- GRanges(
        seqnames = paste0("chr", bed$chrom), 
        ranges = IRanges(bed$chromStart, bed$chromEnd)
    )

    # get overlap
    hits <- findOverlaps(anno_gr, gr)
    arche_cpgs <- ann450k[queryHits(hits), ] |> rownames()
    return(arche_cpgs)
}

get_arche_cpgs <- function(ann450k, anno_gr) {

    # get cpgs in arche regions
    arche1_cpgs <- anno_cpgs("ARCHE1")
    arche2_cpgs <- anno_cpgs("ARCHE2")
    arche3_cpgs <- anno_cpgs("ARCHE3")
    arche4_cpgs <- anno_cpgs("ARCHE4")
    arche5_cpgs <- anno_cpgs("ARCHE5")
    arche6_cpgs <- anno_cpgs("ARCHE6")

    # compile df of cpgs in arches
    arche_cpgs <- data.frame(
        CpGs = c(arche1_cpgs, arche2_cpgs, arche3_cpgs, arche4_cpgs, arche5_cpgs, arche6_cpgs),
        ARCHE = c(
            rep("ARCHE1", length(arche1_cpgs)),
            rep("ARCHE2", length(arche2_cpgs)),
            rep("ARCHE3", length(arche3_cpgs)),
            rep("ARCHE4", length(arche4_cpgs)),
            rep("ARCHE5", length(arche5_cpgs)),
            rep("ARCHE6", length(arche6_cpgs))
        )
    )
    return(arche_cpgs)
}
