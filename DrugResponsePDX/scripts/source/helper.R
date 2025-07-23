#' Standardize PDX Naming
#'
#' Reference defined mapping to standardize PDX names.
#' @param models string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 

# define mapping
mapping <- c("X108099P1" = "108099", 
             "X16720" = "16720",
             "X48602" = "48602",
             "X53782" = "53782",
             "BPTO95" = "BPTO.95",
             "INS_B014" = "INSB014",
             "INS_B019" = "INSB019",
             "REFS032" = "REF032",
             "REFS034P1" = "REF034", 
             "REFS038" = "REF038",
             "REFB047" = "REF047",
             "REF047_S26_L002" = "REF047",
             "REFS036P1A" = "REF036",
             "REF_S_038" = "REFS038"
)

map_pdx <- function(models) {

    for (i in 1:length(models)) {
        model = models[i]
        if (model %in% names(mapping)) {
            models[i] <- unname(mapping[model])
        }
    }
    return(models)
}


#' Extract PDX ARCHE signature scores
#'
#' Append columns to treatment response matrices with PDX ARCHE signature scores.
#' @param df string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 
get_ARCHE <- function(df) {
    for (i in 1:nrow(df)) {
        sample = df$patient.id[i]
        df$ARCHE1[i] <- scores[rownames(scores) == sample,]$Signature1
        df$ARCHE2[i] <- scores[rownames(scores) == sample,]$Signature2
        df$ARCHE3[i] <- scores[rownames(scores) == sample,]$Signature3
        df$ARCHE4[i] <- scores[rownames(scores) == sample,]$Signature4
        df$ARCHE5[i] <- scores[rownames(scores) == sample,]$Signature5
        df$ARCHE6[i] <- scores[rownames(scores) == sample,]$Signature6
    }
    return(df)
}

#' Plot waterfall plot of ARCHE score coloured by mRECIST 
#'
#' Append columns to treatment response matrices with PDX ARCHE signature scores.
#' @param df dataframe. Treatment response dataframe with ARCHE scores concatenated
#' @param arche string. ARCHE name, one of ARCHE[1:6] 
#' @param df string. Drug name (should match entry under 'drug' column of df)
#' @param plot.indiv boolean. TRUE if waterfall plot should be saved as an individual file.
#' @return ggplot object of waterfall plot.
#' 
waterfall_mRECIST <- function(df, arche, drug, plot.indiv = F) {
 
    # create ranking order
    df <- df[order(df[[arche]], decreasing = T),]
    df$rank <- 1:nrow(df)

    # obtain min and max values
    y <- max(abs(min(df[[arche]])), max(df[[arche]])) |> ceiling()

    # create waterfall plot coloured by mRESCIST
    p <- ggplot(df, aes(x = rank, y = .data[[arche]], fill = mRECIST)) + 
        geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
        scale_fill_manual(values = c("#136F63", "#9DCBBA", "#FFADA1", "#B02E0C"),
                        labels = c('CR' = 'Complete\nResponse', 
                                'PD' = 'Progressive\nDisease', 
                                'SD' = 'Stable\nDisease',
                                'PR' = 'Partial\nResponse')) +
        theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylim(c(-y, y)) + labs(x = "PDX Model", y = paste(arche, "Score"), fill = "mRECIST") 
    
    # print out waterfall plot if needed
    if (plot.indiv == TRUE) {
        path = paste0("DrugResponsePDX/results/figures/indiv_plots/mRECIST_waterfall_", arche, "_", drug, ".png")
        png(path, width=175, height=125, units='mm', res = 600, pointsize=80)
        print(p)
        dev.off()
    }
    return(p)
}