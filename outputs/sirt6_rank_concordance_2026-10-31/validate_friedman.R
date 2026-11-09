library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "outputs/sirt6_rank_concordance_2026-10-31"
cases <- fromJSON(file.path(out, "validation_cases.json"), simplifyVector = FALSE)
rows <- lapply(cases, function(z) {
  ranks <- do.call(rbind, lapply(z$ranks, unlist))
  fit <- friedman.test(ranks)
  statistic <- unname(fit$statistic)
  w <- statistic / (nrow(ranks) * (ncol(ranks) - 1))
  stopifnot(abs(statistic - z$expected$friedman_chi2) < 1e-10)
  stopifnot(abs(w - z$expected$kendall_w) < 1e-12)
  stopifnot(abs(log(fit$p.value) - log(z$expected$p_chi2)) < 1e-10)
  data.frame(dataset = z$dataset, group = z$group,
             n_species = nrow(ranks), n_genes = ncol(ranks),
             R_friedman_chi2 = statistic, R_p = fit$p.value, R_kendall_w = w,
             chi2_abs_difference = abs(statistic - z$expected$friedman_chi2),
             W_abs_difference = abs(w - z$expected$kendall_w))
})
validation <- do.call(rbind, rows)
write.table(validation, file.path(out, "R_friedman_validation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(nrow(validation), "tests independently reproduced using R stats::friedman.test.\n")
cat("Maximum statistic difference:", max(validation$chi2_abs_difference), "\n")
cat("Maximum W difference:", max(validation$W_abs_difference), "\n")
