# helper to calculate HR
calculate_HR <- function(df, group, ref) {
    df$ref <- ifelse(df[[ref]] == group, group, "Other")
    df$ref <- relevel(factor(df$ref), ref = "Other")

    fit <- coxph(Surv(time_cutoff, event_cutoff) ~ ref, data = df)
    s <- summary(fit)
    res <- data.frame(
        Group = group,
        HR = s$conf.int[, "exp(coef)"],
        lower = s$conf.int[, "lower .95"],
        upper = s$conf.int[, "upper .95"],
        pvalue = s$coefficients[, "Pr(>|z|)"]
    )
    return(res)
}

# helper function to plot survival plot
indiv_plot_survival <- function(df, group, ref, cutoff, pal) {

    df$ref <- ifelse(df[[ref]] == group, group, "Other")
    df$ref <- factor(df$ref, levels = c("Other", group))
    fit <- survfit(Surv(time_cutoff, event_cutoff) ~ ref, data = df)

    p <- ggsurvplot(
        fit,
        data = df,
        size = 2,
        pval = TRUE,
        risk.table = TRUE,
        palette = c("#979797", pal[[group]]),
        legend.title = "",
        xlab = "Time (days)",
        xlim = c(0, cutoff),
        break.time.by = 365
    )
    p$plot <- p$plot + 
        theme(
            legend.key.size = unit(3, "cm"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            
        )
    p$table <- p$table + theme(plot.title = element_text(size = 12))

    filename <- paste0("data/results/figures/1-Signatures/ARCHEG_survivalplots/", group, "_", as.character(cutoff/365), ".png")
    png(filename, width=4.5, height=5, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()
}

# function to plot survival plots
plot_survival <- function(archeG, pheno, cutoff) {

    cutoff <- cutoff*365
    pheno$time_cutoff <- pmin(pheno$overall_survival, cutoff)
    pheno$event_cutoff <- ifelse(pheno$overall_survival <= cutoff, pheno$status, 0)

    # add to plot
    archeG$time_cutoff <- pheno$time_cutoff[match(archeG$Sample, rownames(pheno))]
    archeG$event_cutoff <- pheno$event_cutoff[match(archeG$Sample, rownames(pheno))]

    HR_res <- data.frame(matrix(nrow=0, ncol=0))

    for (arche in paste0("ARCHE", 1:6)) {
        indiv_plot_survival(archeG, arche, "ARCHEG", cutoff, pal = ARCHE_pal)
        HR_res <- rbind(HR_res, calculate_HR(archeG, arche, "ARCHEG"))
    }

    for (subtype in unique(archeG$PAM50)) {
        indiv_plot_survival(archeG, subtype, "PAM50", cutoff, pal = subtype_pal)
        HR_res <- rbind(HR_res, calculate_HR(archeG, subtype, "PAM50"))
    }
    
        toPlot <- HR_res %>%
            mutate(
                label_hr = sprintf("%.2f (%.2f-%.2f)", HR, lower, upper),
                label_p  = case_when(
                pvalue < 0.001 ~ "p < 0.001",
                TRUE      ~ sprintf("p = %.3f", pvalue)
                ),
                Group = factor(Group, levels = c(rev(names(subtype_pal)[1:5]), paste0("ARCHE", 6:1)))  # keeps order top-to-bottom
            )

    count <- rbind(as.matrix(table(archeG$ARCHEG)), as.matrix(table(archeG$PAM50))) |> as.data.frame()
    colnames(count) <- "N"
    toPlot$n <- count$N[match(toPlot$Group, rownames(count))]


    p1 <- ggplot(toPlot, aes(x = HR, y = Group)) +
        annotate("rect", xmin = 0.3, xmax = 5.2, ymin =2.5, ymax = 3.5, fill = "#CAE2BC", alpha = 0.5) +
        annotate("rect", xmin = 0.3, xmax = 5.2, ymin =5.5, ymax = 6.5, fill = "#CAE2BC", alpha = 0.5) +
        annotate("rect", xmin = 0.3, xmax = 5.2, ymin =6.5, ymax = 7.5, fill = "#C17F7F", alpha = 0.5) +
        annotate("rect", xmin = 0.3, xmax = 5.2, ymin =8.5, ymax = 9.5, fill = "#C17F7F", alpha = 0.5) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = lower, xmax = upper), width = 0.15, linewidth = 0.8) +
        geom_point(aes(size = n, fill = Group), shape = 22, color = "black") +
        scale_size_continuous(range = c(3, 8), guide = "none") +
        scale_fill_manual(values = c(ARCHE_pal, subtype_pal)) +
        scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
        labs(
            x = "Hazard Ratio (95% CI)",
            y = NULL
        ) +
        theme_classic(base_size = 13) +
        theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 12),
            legend.position = "none",
            plot.margin = margin(5, 0, 5, 0)
        )

    p2 <- ggplot(toPlot, aes(x = 0, y = Group, label = label_hr)) +
        annotate("rect", xmin = 0, xmax = 1.5, ymin =2.5, ymax = 3.5, fill = "#CAE2BC", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = 1.5, ymin =5.5, ymax = 6.5, fill = "#CAE2BC", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = 1.5, ymin =6.5, ymax = 7.5, fill = "#C17F7F", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = 1.5, ymin =8.5, ymax = 9.5, fill = "#C17F7F", alpha = 0.5) +
        geom_text(hjust = 0, size = 4.5) + scale_x_continuous(limits = c(0, 1.5), expand = c(0, 0)) + theme_void() 
    p3 <- ggplot(toPlot, aes(x = 0, y = Group, label = label_p)) +
        annotate("rect", xmin = 0, xmax = 0.6, ymin =2.5, ymax = 3.5, fill = "#CAE2BC", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = 0.6, ymin =5.5, ymax = 6.5, fill = "#CAE2BC", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = 0.6, ymin =6.5, ymax = 7.5, fill = "#C17F7F", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = 0.6, ymin =8.5, ymax = 9.5, fill = "#C17F7F", alpha = 0.5) +
        geom_text(hjust = 0, size = 4.5) + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_void() +
        theme(plot.margin = margin(5, 0, 5, 5))

    p <- p1 + p2 + p3 + plot_layout(width = c(2, 1, 1))

    filename <- paste0("data/results/figures/1-Signatures/ARCHEG_survivalplots/HR_", as.character(cutoff/365), ".png")
    ggsave(filename, p, width = 6, height = 4)

}