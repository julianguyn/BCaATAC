#' Get ARCHE assignments from NMF output
#' 
#' Helper function to speed up loading in ARCHE matrix
#' 
get_ARCHE <- function() {

  atac <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)
  atac$ARCHE <- paste0("ARCHE", 1:6)
  mat <- reshape2::melt(atac)
  mat$variable <- gsub("X", "", mat$variable)

  # save signature assignment
  mat$signature_assign <- ""
  for (sample in mat$variable) {
    tmp <- mat[mat$variable == sample,]
    mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$ARCHE)
  }
  mat <- mat[,colnames(mat) %in% c("variable", "signature_assign")] |> unique()

  mat$rank <- 1:nrow(mat)
  return(mat)
}

#' Get ARCHE assignments from NMF output (keep as mat)
#' 
#' Helper function to speed up loading in ARCHE matrix
#' Must have loaded in metadata first
#' 
get_ARCHE_mat <- function() {

  atac <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)
  rownames(atac) <- paste0("ARCHE", 1:6)
  colnames(atac) <- meta$Sample.Name[match(gsub("X", "", colnames(atac)), meta$ATAC.Seq.File.Name)]

  return(atac)
}

#' Get tumour metadata for mutations
#' 
get_meta_mut <- function() {
  meta <- read.csv("MolecularSigAnalysis/data/TCGA_sourcefiles.csv")
  meta$SNV.File.Name <- paste0("MolecularSigAnalysis/data/maf/", meta$SNV.File.Name)
  meta <- meta[-which(meta$SNV.File.Name == "MolecularSigAnalysis/data/maf/NA"),]
  meta$Signature <- gsub("Signature", "ARCHE", meta$Signature)
  return(meta)
}

#' BCa-relevant mutations
#' 
#' List of mutation of interest
#' 
bca_mutations <- c("PIK3CA", "TP53", "CDH1", "ERBB2", "MUC16", "TG", "TTN", "MYC",
                   "QSER1", "SPTA1", "PTPRB", "GATA3", "USH2A", "TMCC3", "NCOR1", 
                   "ZFHX4", "PTEN", "BRCA1", "BRCA2", "ATM", "CHEK2")

#' Run DEG
#' 
#' @param counts dataframe. Counts matrix
#' @param meta dataframe. Metadata table
#' @param arche string. ARCHE to run by
#' @param subtype string. Subtype to run by 
#' 
run_DEG <- function(counts, meta, arche = NA, subtype = NA) {

  # create colData
  colData <- data.frame(
    SampleID = meta$Sample.Name,
    ARCHE = ifelse(meta$Signature == arche, arche, "Other"),
    Subtype = ifelse(meta$Subtype == subtype, subtype, "Other")
  )
  
  if (!is.na(arche)) {  # if running ARCHEx vs all
    colData$ARCHE <- factor(colData$ARCHE, levels = c("Other", arche))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ ARCHE) |>
                                suppressMessages()
    label <- paste0("ARCHE_", arche, "_vs_Other")
  } else {            # if running Subtypex vs all
    colData$Subtype <- factor(colData$Subtype, levels = c("Other", subtype))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ Subtype) |>
                                suppressMessages()
    label <- paste0("Subtype_", subtype, "_vs_Other")
  }

  # DEG
  dds <- DESeq(dds) |> suppressMessages()
  res <- results(dds, name=label)

  write.csv(res, file = paste0("MolecularSigAnalysis/results/data/DEG/", label, ".csv"), quote = F, row.names = T)

  return(res)
}

#' Run single sample gsea
#' 
#' @param counts data.frame. Gene counts matrix
#' @param gmt GeneSetCollection or list. 
#' @param label string. Signature label for saving
#' 
run_ssgsea <- function(counts, gmt, label) {

  es <- gsva(ssgseaParam(counts, gmt), verbose = FALSE)
  write.csv(es, file = paste0("MolecularSigAnalysis/results/data/", label, "_ES.csv"), quote = F, row.names = T)

  return(es)
}

#' Format mutational sig dataframes for plotting
#' 
#' Helper function add ATAC signature to mutational sig dataframe
#' Used in plot_heatmap_mutsig()
#' 
format_sig <- function(cosm_df) {
  cosm_df <- as.data.frame(cosm_df)
  cosm_df$Sample.Name <- rownames(cosm_df)
  cosm_df$Signature <- meta$Signature[match(cosm_df$Sample.Name, meta$Sample.Name)]
  cosm_df$rank <- meta$rank[match(cosm_df$Sample.Name, meta$Sample.Name)]
  cosm_df$Subtype <- meta$Subtype[match(cosm_df$Sample.Name, meta$Sample.Name)]
  cosm_df <- reshape2::melt(cosm_df, id = c("Sample.Name", "Signature", "Subtype", "rank"))
  cosm_df$rank <- factor(cosm_df$rank, levels = max(cosm_df$rank):1)
  return(cosm_df)
}

#' Correlate each pair of mut sig to ARCHE
#' 
#' @param sig data.frame. Output of compareSignatures |> t() or ssgsea
#' @param label string. Label for file output
#' 
corr_signatures <- function(sig, label) {

  # load in ARCHE matrix
  atac <- get_ARCHE_mat()

  # keep common samples
  common <- intersect(colnames(atac), colnames(sig))
  atac <- atac[,match(common, colnames(atac))]
  sig <- sig[,match(common, colnames(sig))]
  
  # initiate dataframe to hold results
  corr <- data.frame(matrix(nrow=0, ncol=3))
  colnames(corr) <- c("ATAC.Sig", "Mol.Sig", "Corr")
  
  # loop through and correlate
  for (i in seq_along(rownames(atac))) {
    
    atac.sig <- rownames(atac)[i]
    
    for (j in seq_along(rownames(sig))) {
      
      molsig <- rownames(sig)[j]
      s <- cor(as.numeric(atac[i,]), as.numeric(sig[j,]), method = "spearman")
      corr <- rbind(corr, data.frame(ATAC.Sig = atac.sig, Mol.Sig = molsig, Corr = s))
    }
  }

  # format for plotting
  if (label == "og30") {
    corr$Mol.Sig <- factor(corr$Mol.Sig, levels = paste0("COSMIC_", 1:30))
  } else if (label == "v3") {
    corr$Mol.Sig <- factor(corr$Mol.Sig, levels = paste0("SBS", 
                           c(1:6, "7a", "7b", "7c", "7d", 8:9, 
                           "10a", "10b", 11:16, "17a", "17b",
                           18:60, 84:85)))
  } else {
    corr$Mol.Sig <- gsub("HALLMARK_", "", corr$Mol.Sig)
  }

  write.csv(corr, file = paste0("MolecularSigAnalysis/results/data/", label, "_corr_atac.csv"), quote = F, row.names = F)
  return(corr)
}
