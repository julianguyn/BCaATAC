#' Function to score ARCHEs from Griffin outputs
#' 
#' @param dir filepath. Directory of Griffin outputs
#' @param sites string. List of site names used on Griffin
#' 
score_arche_cfDNA <- function(dir, sites = NA) {

    # get all griffin files
    files <- list.files(
        dir,
        recursive = TRUE,
        pattern = "GC_corrected.coverage.tsv"
    )
    samples <- gsub("/.*", "", files)

    if (is.na(sites)) sites <- paste0("sig", 1:6, "_top10k")

    # create dataframe to store results
    res <- data.frame(matrix(nrow=0, ncol=3)) 
    colnames(res) <- c("Sample", "ARCHE", "Score")

    for (file in files) {
        df <- fread(paste0(dir, "/", file))
        df <- df[df$site_name %in% sites,]
        df  <- data.frame(Sample = df$sample,
                          ARCHE = df$site_name,
                          Score = 1-df$central_coverage)
        df$ARCHE <- map_griffin(df$ARCHE)
        res <- rbind(res, df)
    }
    return(res)
}
