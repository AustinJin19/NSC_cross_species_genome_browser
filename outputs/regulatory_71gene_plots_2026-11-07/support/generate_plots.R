#!/usr/bin/env Rscript
# Recreate all plots offline from the saved inputs in this support directory.
arg <- grep('^--file=', commandArgs(), value = TRUE)
if (!length(arg)) stop('Run this script using Rscript.')
here <- dirname(normalizePath(sub('^--file=', '', arg[1])))
out <- dirname(here)
source(file.path(here, 'plot_regulatory_mls.R'))
data <- jsonlite::fromJSON(file.path(here, 'species_medians.json'))
genes <- readLines(file.path(here, 'genes.txt'))
tests <- read.delim(file.path(here, 'peaks500kb_without_body_mass_all_tested.tsv.gz'))
mass_tests <- read.delim(file.path(here, 'peaks500kb_with_body_mass_all_tested.tsv.gz'))
results <- vector('list', length(genes))
exclusions <- vector('list', length(genes))
fits <- vector('list', length(genes))
for (i in seq_along(genes)) {
  gene <- genes[i]
  d <- data[data$gene == gene, , drop = FALSE]
  filename <- paste0(gene, '_500kb_regulatory_MLS.png')
  fit <- plot_regulatory_mls(d, file.path(out, filename),
    fdr_results = tests, body_mass_fdr_results = mass_tests, width = 6, height = 6)
  u <- fit$statistics
  m <- fit$body_mass_corrected$statistics
  # Match the exact species sets as well as P values to the full-family scan.
  stopifnot(identical(sort(fit$data$species),
    sort(strsplit(tests$species[match(gene, tests$gene)], ';', fixed = TRUE)[[1]])),
    identical(sort(fit$body_mass_corrected$data$species),
    sort(strsplit(mass_tests$species[match(gene, mass_tests$gene)], ';', fixed = TRUE)[[1]])))
  results[[i]] <- data.frame(gene=gene, plot=filename,
    n_species=u$n, Pearson_r=u$r, nominal_P=u$p_value, BH_adjusted_P=u$BH_FDR,
    n_species_body_mass_corrected=m$n, partial_r_body_mass_corrected=m$r,
    nominal_P_body_mass_corrected=m$p_value, BH_adjusted_P_body_mass_corrected=m$BH_FDR,
    genes_tested=u$FDR_genes_tested, genes_tested_body_mass_corrected=m$FDR_genes_tested,
    included_species=paste(fit$data$species, collapse=';'),
    included_species_body_mass_corrected=paste(fit$body_mass_corrected$data$species, collapse=';'))
  exclusions[[i]] <- fit$excluded
  fits[[i]] <- list(gene=gene, coefficients=stats::coef(fit$model), statistics=u,
    body_mass_corrected_coefficients=stats::coef(fit$body_mass_corrected$model),
    body_mass_corrected_statistics=m)
  if (i %% 10 == 0 || i == length(genes)) cat(i, '/', length(genes), 'plots generated\n')
}
summary <- do.call(rbind, results)
write.csv(summary, file.path(out, 'plot_statistics.csv'), row.names=FALSE, na='')
write.csv(do.call(rbind, exclusions), file.path(here, 'excluded_gene_species.csv'), row.names=FALSE, na='')
saveRDS(fits, file.path(here, 'fitted_models.rds'))
stopifnot(nrow(summary)==71, all(file.exists(file.path(out, summary$plot))))
cat('BH-adjusted P < 0.05:', sum(summary$BH_adjusted_P < .05), 'uncorrected;',
  sum(summary$BH_adjusted_P_body_mass_corrected < .05), 'body-mass corrected\n')
