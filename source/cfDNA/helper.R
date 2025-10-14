#' Compute ARCHE scores
#'
#' Computes ARCHE score (1 - central coverage) from Griffin output.
#' @param files string. List of Griffin output files.
#' @param label string. Group name (baseline/control/progression)
#' @return A dataframe of ARCHE scores.
#' 
score_ARCHE <- function(files, label) {

    # create dataframe to store results
    res <- data.frame(matrix(nrow=0, ncol=3)) 
    colnames(res) <- c("Sample", "ARCHE", "Score")

    # get baseline results (30)
    for (file in files) {
        df <- fread(file)

        # specify sites to keep
        keep <- c("sig1_top10k", "sig2_top10k", "sig3_top10k", "sig4_top10k", "sig5_top10k", "sig6_top10k")
        df <- df[df$site_name %in% keep,]

        # extract information from griffin output
        df  <- data.frame(Sample = df$sample,
                          ARCHE = paste0("ARCHE", 1:6),
                          Score = 1-df$central_coverage)
        res <- rbind(res, df)
    } 
    res$label <- label
    return(res)
}

#' Get ARCHE summary stats
#'
#' Compute mean and SE of ARCHE scores
#' @param df dataframe. ARCHE scores from score_ARCHE()
#' @param TF_group string. Dataframe TF_group label
#' @return A dataframe of ARCHE score summary stats.
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