#' Function to score ARCHEs from Griffin outputs
#' 
#' @param dir filepath. Directory of Griffin outputs
#' @param meta metadata file
#' 
score_arche_cfDNA <- function(dir, meta = NULL) {

    dir <- paste0("data/rawdata/cfDNA/", dir)

    # get all griffin files
    files <- list.files(
        dir,
        recursive = TRUE,
        pattern = "GC_corrected.coverage.tsv"
    )
    samples <- gsub("/.*", "", files)

    # create dataframe to store results
    res <- data.frame(matrix(nrow=0, ncol=3)) 
    colnames(res) <- c("Sample", "ARCHE", "Score")

    for (file in files) {
        df <- fread(paste0(dir, "/", file))
        df  <- data.frame(Sample = df$sample,
                          Label = df$site_name,
                          Score = 1-df$central_coverage)
        res <- rbind(res, df)
    }

    res$ARCHE <- sub("_.*", "", res$Label)
    res$Subset <- sub(".*_", "", res$Label)

    # add metadata variables
    if (dir != "data/rawdata/cfDNA/preclinical-ARCHE") {
        res$label <- meta$time_id[match(res$Sample, meta$sample_id)]
        res$TF <- meta$metrics_tf[match(res$Sample, meta$sample_id)]
        res$pheno <- meta$pheno_id[match(res$Sample, meta$sample_id)]
        res$TF_group <- ifelse(res$TF > 10, "TF>10", "TF<10")
    }

    return(res)
}

#' Function to summarize ARCHE scores
#' 
#' Compute mean and SE of ARCHE scores
#' @param scores df. Output of score_arche_cfDNA()
#' 
summary_ARCHE <- function(df, TF_group) {
    summary_df <- df %>%
        group_by(ARCHE, label) %>%
        summarise(
            mean_score = mean(Score, na.rm = TRUE),
            se_score = sd(Score, na.rm = TRUE) / sqrt(n())
        ) %>%
        ungroup()
    summary_df <- as.data.frame(summary_df)
    summary_df$TF_group <- TF_group
    return(summary_df)
}
