#' Correlate each pair of molecular sig to ARCHE
#' 
#' @param sig data.frame. Output of compareSignatures |> t() or ssgsea
#' @param label string. Label for file output
#' @param data string. "tumour" or "ccls"
#' 
corr_signatures <- function(sig, label, data = "tumour") {

  # load in ARCHE matrix
  if (data == "tumour") {
    atac <- get_arche_tcga(mat=TRUE)
  } else if (data == "ccls") {
    atac <- get_arche_cells()
  }

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

  write.csv(corr, file = paste0("data/results/data/2-MolecularSigAnalysis/molecularsig/", label, "_corr_atac.csv"), quote = F, row.names = F)
  return(corr)
}
