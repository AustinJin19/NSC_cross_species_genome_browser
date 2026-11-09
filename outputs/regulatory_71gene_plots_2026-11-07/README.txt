71 regulatory activity / MLS plots

One PNG per gene from the original 71 positive TSS +/-500 kb regulatory hits
in Regulatory_500kb_MLS_and_GO_2026-10-19.xlsx. The gene list is unchanged.

Each image is square: 6 x 6 inches, 1800 x 1800 pixels at 300 dpi.
X: log10(maximum lifespan in years).
Y: log10(species-median peak count within TSS +/-500 kb).
Black points show species medians without species-name labels.
The red dashed line, r and nominal P describe the uncorrected ordinary fit.
Both BH-adjusted P values are shown, with and without log10(body mass) as a
covariate. BH is calculated separately across all 18,621 tested genes in each
current model, NOT across these selected 71 genes.

The data and statistical method match the current SIRT6 plot dated 2026-11-07.
Human, macaque and bowhead whale are excluded. MLS lower bounds are excluded.
The highest species-median locus is chosen before eligibility filtering.
Species with <50% of the requested genomic window available are excluded.
No positive-ATAC-signal coverage threshold is used. Log transforms have no
pseudocount. Missing/nonpositive body mass is omitted from the corrected model.
This is ordinary regression, not PGLS.

plot_statistics.csv lists each plot and its exact P values, BH-adjusted P
values, correlation coefficients, sample sizes and included species.
All 71 original genes are plotted. With current data/name mapping and BH
families, 51 have uncorrected BH P < 0.05 and 70 have body-mass-corrected
BH P < 0.05. DENND10 has corrected BH P = 0.0555850.

The support folder contains the selected species measurements, both complete
BH testing families, model coefficients, provenance and scripts to reproduce the
plots. From any directory, use:

  Rscript /full/path/to/this/folder/support/generate_plots.R

R packages required: ggplot2 and jsonlite. The plotting script runs offline.
prepare_inputs.py reconstructs inputs from the validated current project
cache and refuses to combine changed source files with the saved BH families.
