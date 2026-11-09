#!/usr/bin/env Rscript
# Run from any working directory: Rscript /full/path/to/run_example.R
# Reproduces the example from the saved 2026-11-07 browser measurements,
# without needing a running website. Uses ggplot2.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("Run this example with Rscript.")
here <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
project <- normalizePath(file.path(here, "../.."))
source(file.path(project, "cross_species_visual_tool/plot_regulatory_mls.R"))
d <- readRDS(file.path(here, "SIRT6_500kb_species_medians.rds"))
tests_unadjusted <- read.delim(file.path(here, "peaks500kb_without_body_mass_all_tested.tsv.gz"))
tests_adjusted <- read.delim(file.path(here, "peaks500kb_with_body_mass_all_tested.tsv.gz"))
unadjusted <- plot_regulatory_mls(d, file.path(here, "SIRT6_500kb_regulatory_MLS.png"),
  fdr_results = tests_unadjusted, body_mass_fdr_results = tests_adjusted)
adjusted <- plot_regulatory_mls(d,
  file.path(here, "SIRT6_500kb_regulatory_MLS_mass_corrected.png"),
  body_mass_correction = TRUE, fdr_results = tests_adjusted)
print(unadjusted$statistics)
print(adjusted$statistics)
# For vector PDF instead, change the output extension to .pdf.
# Use show_statistics = FALSE for axes, points and fitted line only.
