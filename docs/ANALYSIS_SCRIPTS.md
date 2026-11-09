# Analysis scripts and regulatory plots

The repository preserves analysis source under its original project-relative
paths, including `output/`, `outputs/`, and `tmm_pca/`. These directories contain
code in Git, not the local result datasets. Historical analyses record the
settings used at the time and may differ from the current website. Some scripts
also use local absolute paths or the document-generation runtime in which they
were written; check their inputs before running on another machine.

The repository includes compact browser metadata: the gene-name mapping,
additional species configuration, lifespan/body-mass values, sample-QC policy,
TimeTree, and its phylogenetic covariance. ATAC tracks, peaks, transcript
indexes, SQLite databases, reference genomes, analysis caches, figures and
workbooks must be transferred separately. See [data requirements](DATA_REQUIREMENTS.md).

## Regulatory activity plots

Source `cross_species_visual_tool/plot_regulatory_mls.R` in R. It requires
`ggplot2`; loading browser measurements also requires `jsonlite`.

```r
source("cross_species_visual_tool/plot_regulatory_mls.R")
d <- load_regulatory_activity("SIRT6", window_kb = 500, measurement = "count")
tests <- read.delim("outputs/regulatory_mls_plot_2026-11-07/peaks500kb_without_body_mass_all_tested.tsv.gz")
mass_tests <- read.delim("outputs/regulatory_mls_plot_2026-11-07/peaks500kb_with_body_mass_all_tested.tsv.gz")
plot_regulatory_mls(d, "SIRT6_500kb.png", fdr_results = tests,
                    body_mass_fdr_results = mass_tests)
```

The default export is 6 by 6 inches, with black points and a red dashed
uncorrected fitted line. It shows log10 median peak count versus log10 MLS,
the uncorrected nominal P, and separate BH-adjusted P values from models with
and without log10 body mass. Human, macaque and bowhead whale are excluded by
default. Genomic window availability must be at least 50%; no positive-signal
coverage threshold is applied. This function uses ordinary regression.

The two complete all-tested tables must correspond to the same measurements,
window, species policy and model. The function checks the plotted gene's P
and sample size against each family and recomputes BH across **all tested
genes**, not just the selected 71. Those tables are generated data and are not
stored in Git.

`outputs/regulatory_mls_plot_2026-11-07/recalculate_fdr.py` creates the current
500-kb testing families from installed project data, using the vectorized
counter in `outputs/regulatory_windows_2026-10-22/analyze_windows.py`.
`run_example.R` reproduces the saved SIRT6 example when its RDS and testing
tables have been restored to that folder.

## Original 71-gene batch

The batch source and gene list are in
`outputs/regulatory_71gene_plots_2026-11-07/support/`.

- `prepare_inputs.py` extracts the selected species medians from the validated
  quantification cache. It needs the original 71-gene list and `all_tested.tsv`
  from `outputs/regulatory-500kb-mass-2026-10-19/`, the current browser metadata,
  and the current 500-kb scan outputs. Changed source fingerprints require a
  new scan before plotting.
- `generate_plots.R` runs offline after the prepared species-measurement JSON
  and both full testing tables have been restored/generated in `support/`.
  It exports one square PNG per gene and `plot_statistics.csv` with both models.
- The local November 7 batch retains all original 71 genes. With that analysis
  version, 70 meet body-mass-corrected BH P < 0.05; DENND10 has BH P 0.055585.

## ATAC super-enhancers

`pipelines/atac_super_enhancers/` contains the filtered-BAM workflow, a manifest
template, reference sequence dictionaries, and pinned ROSE source with its
upstream license and checksums. BAM locations must be supplied on the cluster.
Follow its [README](../pipelines/atac_super_enhancers/README.md) for environment,
validation, one-sample testing, SLURM submission and per-species summaries.
