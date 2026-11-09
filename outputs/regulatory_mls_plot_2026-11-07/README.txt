Standalone R function: ../../cross_species_visual_tool/plot_regulatory_mls.R
Black points, dashed red fitted line, no species labels. White background, light grid,
black panel border. Gene name is the top title, followed by nominal P and BH FDR.
Export formats are selected by filename (e.g. PNG or PDF).

Dependencies: ggplot2; jsonlite for loading new data from the running browser.
Source the function file, then:

  d <- load_regulatory_activity("SIRT6", window_kb = 500, measurement = "count")
  tests <- read.delim("peaks500kb_without_body_mass_all_tested.tsv.gz")
  result <- plot_regulatory_mls(d, "SIRT6_regulatory_MLS.pdf", fdr_results = tests)

To show both uncorrected and body-mass-corrected BH-adjusted P values on the
same uncorrected scatter plot:

  mass_tests <- read.delim("peaks500kb_with_body_mass_all_tested.tsv.gz")
  result <- plot_regulatory_mls(d, "SIRT6_both_FDR.png", fdr_results = tests,
    body_mass_fdr_results = mass_tests)

The dashed red line, r and nominal P describe the uncorrected fit. The additional
BH-adjusted P comes from a separately fitted model with log10(body mass) as a
covariate, checked against its own full testing family. Its fitted model,
included/excluded species and statistics are returned as result$body_mass_corrected.
If its species count differs, that count is also printed on the plot.
The default SIRT6 example now shows both BH values: 0.0174 uncorrected and
0.0345 body-mass corrected. The separate corrected example uses the corrected
fitted line and nominal P; it is explicitly labelled as a corrected fit.

Use measurement = "signal" for integrated CPM in peaks. Available windows:
20, 50, 100, 200 or 500 kb on each side of the TSS.
Use body_mass_correction = TRUE for the log10(body mass) covariate, together
with the separate peaks500kb_with_body_mass_all_tested.tsv.gz FDR family.
Use log_activity = FALSE and/or log_mls = FALSE for raw axes.
Use show_statistics = FALSE to omit the numerical annotation.
The gene title is inferred from data$gene; for custom input use gene = "SIRT6".
Default export size is 6 x 6 inches (1:1 canvas; 1800 x 1800 pixels at 300 dpi).
Use width and height (inches) to change export size; dpi defaults to 300.
The function also accepts your own data frame of species, activity (raw median),
and mls (years), plus body_mass_g when mass adjustment is enabled.

Default axes: log10(species-median regulatory activity) vs log10(MLS).
Ordinary least squares; independent species, not PGLS. The plot shows nominal
P and BH-adjusted P across all valid genes for that same window and model.
With mass correction, the fitted line holds log10(body mass) at its included-
species mean. The points remain observed activities on the chosen axes.
Reported adjusted r is the partial correlation, not the simple scatter r.
The returned result contains plot, model, data, excluded, statistics and
fdr_results. Without fdr_results, adjusted P is reported as unavailable.

FDR input must include ALL tested genes, including nonsignificant genes and
both correlation signs, with one row per gene and columns gene and p_value.
NA P values for untestable genes are omitted from the adjustment family.
The function recomputes BH with R's p.adjust and verifies the plotted gene's
nominal P matches the input table. It also checks model and window metadata
when supplied. Do not reuse these +/-500 kb log/log peak-count FDR tables
for another window, raw axes, CPM signal in peaks, PGLS or different filters.

recalculate_fdr.py rebuilds the complete +/-500 kb peak-count testing families
from current browser inputs. It uses the current gene-name resolver and sample
catalog, checks vectorized counts against the production scorer, and requires
at least five eligible species per gene. BH is applied separately to the two
models, across all valid tested genes. See fdr_manifest.json for provenance.

Data loading follows the browser's highest-median-locus selection. For these
plots, exclude human, macaque, bowhead, missing measurements, MLS lower bounds, and
species with <50% of the genomic window available. Exactly 50% is retained.
No positive-signal coverage threshold is applied.
Nonpositive activity is excluded only on the log activity scale. Missing or
nonpositive body mass is excluded only when mass correction is enabled.
For manually supplied data, region/lower-bound filtering is available when
the corresponding columns documented in the function file are provided.

Examples use SIRT6 TSS +/-500 kb peak counts from the local browser on
2026-11-07. The saved RDS contains all 34 species, including exclusion fields.
Both example models include 26 species. Run run_example.R to reproduce PNGs
offline. Excluded: human, macaque, bowhead whale, cave nectar bat (MLS lower bound),
Bramus lutescens, muskrat, eastern mole and fox squirrel (region availability).
Regression slopes and P values were independently checked using Python/SciPy.
