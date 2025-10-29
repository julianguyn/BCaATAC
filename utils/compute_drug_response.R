#' List of functions:
#' format_combinations()
#' compute_ci()
#' compute_pc()
#' compute_meta()

#' FDR correction and format dataframes 
#' 
#' Helper function for compute_ci() and compute_pc()
#' @param combinations dataframe. Intermediary from computeCI() and computePC().
#' @return A dataframe with FDR correction and additional formating.
#' 
format_combinations <- function(combinations, label) {
    
    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    #combinations$FDR <- p.adjust(combinations$pvalue, method = "BH", n = length(combinations$pvalue))
    combinations$FDR <- ave(combinations$pval, combinations$drug, FUN = function(p) p.adjust(p, method = "BH"))
    combinations$FDRsig <- ifelse(combinations$FDR < 0.05, TRUE, FALSE)

    # format dataframe for plotting (ordered by CI/PC already)
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$signature, "_", combinations$drug)
    combinations$pset <- c(rep(label, nrow(combinations)))

    return(combinations)
}

#' Compute concordance index
#'
#' Computes concordance index between ARCHE scores and drug AAC.
#' @param signature_scores dataframe. ARCHE scores from get_scores().
#' @param sensitivity_data dataframe. Drug sensitivity from get_drugsen().
#' @param label string. PSet name
#' @return A dataframe of CI and metrics.
#' 
compute_ci <- function(signature_scores, sensitivity_data, label) {

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("signature", "drug", "ci", "pvalue", "se", "upper", "lower")
    combinations$signature <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){

        ci <- survcomp::concordance.index(as.numeric(sensitivity_data[combinations$drug[i],]), # drug AAC
                        surv.time = as.numeric(unlist(-signature_scores[combinations$signature[i],])), # ARCHE score
                        surv.event = rep(1,length(sensitivity_data)), 
                        outx = TRUE, method="noether", na.rm = TRUE)

        combinations$pvalue[i] <- ci$p.value
        combinations$ci[i] <- ci$c.index
        combinations$se[i] <- ci$se
        combinations$upper[i] <- ci$upper
        combinations$lower[i] <- ci$lower
    }

    # format dataframe and FDR correction
    combinations <- combinations[order(combinations$ci),]
    combinations <- format_combinations(combinations, label)
    
    return(combinations)
}

#' Compute Pearson's correlation
#'
#' Computes PC between ARCHE scores and drug AAC.
#' @param signature_scores dataframe. ARCHE scores from get_scores().
#' @param sensitivity_data dataframe. Drug sensitivity from get_drugsen().
#' @param label string. PSet name
#' @return A dataframe of CI and metrics.
#' 
compute_pc <- function(signature_scores, sensitivity_data, label) {
    
    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("signature", "drug", "pc", "pvalue", "n", "upper", "lower")
    combinations$signature <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){
        AAC <- as.numeric(sensitivity_data[combinations$drug[i],])
        if (length(AAC[!is.na(AAC)]) > 3) {
            pc <- cor.test(AAC, # drug AAC
                    as.numeric(unlist(signature_scores[combinations$signature[i],])), # ARCHE score
                    method = 'pearson', alternative = 'two.sided')

            combinations$pvalue[i] <- pc$p.value
            combinations$pc[i] <- pc$estimate
            combinations$n[i] <- length(AAC)
            combinations$upper[i] <- pc$conf.int[1]
            combinations$lower[i] <- pc$conf.int[2]
        }
    }

    # format dataframe and FDR correction
    combinations <- combinations[order(combinations$pc),]
    combinations <- format_combinations(combinations, label)
    
    return(combinations)
}

#' Perform meta analysis
#'
#' Compute meta estimates for signature drug response associations in >=3 PSets.
#' @param df dataframe. Dataframe from computePC()
#' @return A dataframe of meta analysis results.
#' 
compute_meta <- function(df) {
    
    # keep only signature-drug pairs that are present in at least 3 PSets
    df <- df[which(df$pairs %in% names(table(df$pairs)[table(df$pairs) > 2])),]

    # data frame to hold meta estimates
    estimates <- as.data.frame(matrix(data = NA, nrow = length(unique(df$pairs)), ncol = 10))
    colnames(estimates) <- c("signature","drug", "pair", "TE", "seTE", "upper", "lower", "pval", "pval.Q", "I2")

    # perform meta-analysis
    for (i in 1:length(unique(df$pair))) {
        
        pair <- unique(df$pairs)[i]
        tmp <- df[which(df$pair == pair),]
        meta <- metacor(pc, n, data = tmp, method.tau = "DL", studlab = tmp$pset)
        
        estimates$signature[i] <- tmp$signature[1]
        estimates$drug[i] <- tmp$drug[1]
        estimates$pair[i] <- pair
        estimates$TE[i] <- meta$TE.random
        estimates$seTE[i] <- meta$seTE.random
        estimates$upper[i] <- meta$upper.random
        estimates$lower[i] <- meta$lower.random
        estimates$pval[i] <- meta$pval.random
        estimates$pval.Q[i] <- meta$pval.Q
        estimates$I2[i] <- meta$I2

        if (abs(meta$TE.random) > 0.4) {
            # plot forest plot
            fileName = paste0("data/results/figures/4-DrugResponse/ClassB/meta/",pair,".png")
            
            png(fileName, width = 10, height = 4, res = 600, units = "in")
            title <- pair
            forest(meta,
                leftcols = c("studlab", "TE", "seTE", "lower", "upper", "pval"),
                leftlabs = c(title, "Effect", "SE", "95% CI \n Lower", "95% CI \n Upper", "P value"),
                xlab = "effect estimate", lab.e = "Intervention", sortvar = TE, smlab = " ", text.random = "Random effect", 
                print.I2.ci = FALSE, print.Q = TRUE, print.pval.Q = TRUE, digits.sd = 2, print.I2 = TRUE, print.tau2 = TRUE,
                text.random.w = TRUE, colgap.forest.left = "0.5cm", layout = "RevMan5", test.overall.random = TRUE,
                test.overall.common = TRUE, xlim = "symmetric", col.square = "grey70", col.inside = "grey70", col.square.lines = "grey30", 
                col.diamond.random = "#526863", col.diamond.common  = "#BD6B73", ff.xlab = "bold", fontsize = 11, fs.heading = 11.5,
                squaresize = 0.55, scientific.pval = TRUE, lty.random = NULL, lty.fixed  = NULL)
            dev.off()
        }
    }

    estimates$FDR <- ave(estimates$pval, estimates$drug, FUN = function(p) p.adjust(p, method = "BH"))
    
    return(estimates)
}
