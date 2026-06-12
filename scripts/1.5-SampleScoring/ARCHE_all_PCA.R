# -------- RUN ON H4H
library(SummarizedExperiment)

files <- c(
    "/cluster/projects/bhklab/projects/BCaATAC/peak-set-scoring/data/results/TCGA_cell_PDX/chromvar/cell_pdx_tcga.counts_filtered.rds",
    "/cluster/projects/bhklab/projects/BCaATAC/peak-set-scoring/data/results/TCGA_cell_PDX_PDO/chromvar/cell_pdx_pdo_tcga.counts_filtered.rds"
)

for (file in files) {
    print(file)

    outfile <- paste0("results/", sub("\\..*", "", sub(".*chromvar/", "", file)), ".RData")

    counts <- readRDS(file)
    tt <- as.matrix(assay(counts))
    tt <- as.data.frame(t(tt))

    # run PCA
    print("running PCA")
    pca_res <- prcomp(tt)
    samples <- rownames(tt)
    print(head(samples))

    # get variances
    variance_pca <- pca_res$sdev^2
    variance_pca / sum(variance_pca) -> prop_var
    print(head(prop_var))

    pca_res <- pca_res$x[,c(1:10)]
    print(head(pca_res))

    print(paste("Saving to", outfile))
    save(pca_res, samples, prop_var, file = outfile)
}

print("done")
