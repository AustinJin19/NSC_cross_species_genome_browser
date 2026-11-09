# Gene and regulatory activity versus lifespan

The browser has two independent panels. Both follow the searched gene and can use all ATAC species or the species selected in the sidebar. They do not change when the genome view is panned or zoomed.

## Measurements

**Gene activity & lifespan** integrates source bigWig signal over one of these regions:

- Gene body plus 2 kb upstream, with the promoter placed according to strand.
- TSS ±2 kb, ±5 kb, ±100 kb, ±200 kb or ±500 kb. The default is ±200 kb.

**Regulatory activity** provides two measurements within TSS ±20 kb, ±50 kb, ±100 kb, ±200 kb or ±500 kb:

- **Peak count:** the number of sample BED/narrowPeak intervals overlapping the window. Each peak counts once even when only part of it overlaps.
- **CPM signal in peaks:** the integral of bigWig signal over the union of that sample's peaks intersecting the window. Overlapping, nested and touching peaks are merged, and every base contributes once. Signal outside peaks is excluded. The unit is CPM·bp, not a peak count or mean CPM. There is no distance weighting or division by peak length.

All windows are half-open `[TSS - flank, TSS + flank)`. A ±100 kb window spans 200 kb in total. TSS is the outer boundary of the merged gene locus: body start on the positive strand and body end minus one on the negative strand. This can differ from the transcript chosen for the genome-browser anchor. Windows are clipped at chromosome boundaries.

These are bigWig-based scores, not Signac fragment counts or RNA expression. Human tracks retain source-scaled·bp units and are not relabeled as CPM. Unrepresented bases on an available bigWig chromosome contribute zero. A missing file, missing chromosome or incompatible peak assembly is unavailable, not zero. A valid window containing no peaks has score zero.

## Species summaries and correlations

Each sample is scored separately. The plotted value is the median across available samples within each species. When a symbol has multiple loci, the locus with the highest finite species median for the selected measurement is used; ties are resolved by locus identifier. The per-sample table retains all loci and identifies which was selected. Thus each species contributes at most one point.

### Minimum available region: 50%

For every measurement, available fraction = length of the requested region intersected with its chromosome/scaffold bounds / full requested length. Bounds come from the matching bigWig chromosome header. For gene-body activity, the denominator is annotated gene length plus 2,000 bp; for TSS ±N windows, it is 2N. **A species is excluded from the plot and correlation when this fraction is below 0.5; exactly 0.5 is retained.** The check applies after selecting the highest-activity locus, without substituting a lower-activity locus. It is per gene and method, not a permanent species exclusion.

This fraction measures genomic span, not bigWig read coverage, peak coverage or signal intensity; internal FASTA N gaps are not subtracted. CPM-in-peaks uses the whole requested TSS window for this check, not the subset in peaks. All samples contributing to the median must meet the threshold; the table reports their minimum available fraction if chromosome sizes differ. Missing bounds cannot establish eligibility. Scores are not rescaled to extrapolate beyond scaffold ends.

Both tables retain excluded measurements and show available percentage and fit status. API records preserve raw `median` and `samples`, and add `requested_bp`, `available_bp`, `available_fraction` and `region_eligible`; consumers calculating correlations must require `region_eligible == true` on the selected locus. Per-sample `window_details` also contains these fields. The body and peak-count genome-wide analysis scripts enforce the same filter, including when reusing cached body scores. Existing exports and enrichment results are historical and are not automatically regenerated.

Example: muskrat SIRT6 is on a 40,717-bp scaffold. TSS ±50 kb has 40.717% available and is excluded; ±20 kb has 86.3425% available and is eligible. Its gene-body + 2 kb promoter region is fully available.

Both panels offer raw median activity or log10(median activity), and MLS in years or log10(MLS). Log transformation is applied after the species median; nonpositive values are excluded without a pseudocount. Per-sample table values remain untransformed.

Human and macaque are excluded from correlations, while their values remain in the table. The bowhead checkbox controls its exclusion independently in each panel; regulatory activity excludes it by default. Unavailable MLS and lower-bound-only MLS entries are excluded.

The gene-body method fits:

```text
selected activity scale = intercept + beta_MLS * selected MLS scale
                         + beta_length * log10(annotated gene-body length in kb)
```

Length excludes the 2 kb promoter. By default the model reports the partial correlation controlling for log10(length), its two-sided p-value with n−3 degrees of freedom, and partial R². At least five complete species and a non-degenerate design are required. The plotted points retain the selected activity scale, and the line holds log10(length) at its species mean.

All fixed TSS windows, peak counts and signal-in-peaks measurements default to unadjusted Pearson correlation and an ordinary linear fit. At least three species with variable lifespan and activity are required. Browser p-values are nominal, with no phylogenetic or genome-wide multiple-testing correction.

### Optional body mass correction

Each panel has an independent **Body mass correction** checkbox, off by default. Enabling it adds `log10(adult body mass in grams)` to the model on the selected activity and lifespan scales. The gene-body method retains its log10(length) covariate; other methods adjust for log10(mass) alone. Raw measurements are unchanged. The log mass is calculated directly from the displayed positive body-mass value, without a pseudocount. Species with missing, nonfinite or nonpositive mass are excluded only when the option is enabled; all existing region-availability and species exclusions still apply.

```text
Gene body: activity = intercept + beta_MLS * MLS_axis
                     + beta_length * log10(gene length in kb)
                     + beta_mass * log10(body mass in g)
Other methods: activity = intercept + beta_MLS * MLS_axis
                         + beta_mass * log10(body mass in g)
```

The reported partial correlation removes the selected covariates from both the activity response and the lifespan predictor. The two-sided t test uses `df = n - k - 2`, where k is the number of covariates: n−3 for mass alone and n−4 for length plus mass. Mass-only fits require at least four complete species, and length-plus-mass fits require at least five. Constant covariates, collinear predictors or zero residual variation produce no fit; the code does not silently fall back to an unadjusted correlation. The regression line holds covariates at the included species' mean log values, while points remain on the selected activity scale. Adjusted statistics are labeled partial r and partial R².

Gene-symbol matching does not establish one-to-one orthology. Nearby genes can share peaks and signal. Peak calling, sequencing depth, normalization, annotation quality and phylogeny can affect these exploratory comparisons; these measurements do not establish enhancer-to-gene links.

## Required local data

The Git repository contains code only. Supply these generated files and their matching source data separately:

- `cross_species_visual_tool/web_app/data/gene_activity.sqlite`: loci, sample identifiers and original gene-body scores, plus optional `activity_additions` for new species.
- `cross_species_visual_tool/web_app/data/species_mls.json`: MLS, body mass and source provenance.
- `cross_species_visual_tool/web_app/data/additional_species.json`: optional species catalog generated by `import_liver_bundle.py`, needed to reproduce the expanded catalog.
- Matching bigWigs in the folders configured by `server.SPECS` and the additional catalog; matching portable transcript indexes in `cross_species_visual_tool/index/`.
- BED/narrowPeak files for regulatory activity, normally next to each bigWig. `server.peak_path` resolves filenames. ASM3Liver and ASM4Liver use the corrected peaks in `data/atac/liver_atac_bundle/Acomys_cahirinus/`, matching the active ASM bigWigs.

Copy the generated activity database, catalog and trait JSON files when deploying. Copying code alone does not transfer these data. Cached on-demand scores are held in memory and keyed by source file size/modification time and queried coordinates.

## Build and deploy

Run from `cross_species_visual_tool/web_app` with the project dependencies installed:

```bash
python calculate_gene_activity.py --output ../../output/gene_activity_signac_regions
python activity.py --source ../../output/gene_activity_signac_regions
python activity.py --register-additions
python build_mls.py
python server.py --port 8765
```

The first command computes whole-genome gene-body scores and can take substantial time. Existing validated scores and their database can be copied instead. Rebuilding the base database replaces it, so run registration again afterward if needed. Registration preserves the base rows and indexes additional species for on-demand calculation. Rebuild scores and indexes when source annotations or tracks change.

`build_mls.py` reads the local omics database and supplied AnAge archive. It preserves existing verified values when the optional lab CSV is unavailable. When present, `data/metadata/species_traits_user_2026-10-15.tsv` is applied by `apply_user_traits.py`: exact user-supplied values take precedence, with the requested AnAge fallbacks for ASM and ambiguous/missing entries. Keep provenance in the resulting JSON. Body mass is adult weight in grams.

## Optional genome-wide analyses

Install `requirements-dev.txt` for the analysis scripts, including SciPy. These scripts write local output files, which are not committed to Git.

`analyze_regulatory_strength.py --windows 100 200 500` tests peak counts against raw/log10 MLS with raw/log10 species medians. The window list also accepts 20 and 50. It excludes human, macaque and bowhead, requires at least five species, and applies Benjamini–Hochberg correction separately for each window/scale/axis combination. It does not run the CPM-in-peaks method.

`analyze_body_mls_length_mass.py --log-activity --log-length --exclude-bowhead --no-mass` selects the current body-model transformations and exclusions for offline analysis; inspect its other arguments for input/output paths. Earlier exports using other covariates remain separate. Changes to browser settings do not recalculate existing exports or enrichment results.

For the body model with **Body mass correction** enabled, omit `--no-mass`. The peak-count genome-wide script remains an unadjusted analysis; this browser option does not change its exports.

New filtered runs use output folders ending in `_region50` (peak counts: `output/regulatory_strength_region50/`). Earlier unfiltered exports are preserved. Cached body measurements can be reused, but region eligibility is checked again against current chromosome bounds before each correlation.

## Checks

```bash
python -m unittest test_region_availability test_peak_activity test_activity_index test_gene_activity
node test_activity_panel.js
node test_statistics.js
node test_gene_layout.js
```

The peak-activity tests compare interval integration with an independent base-by-base reference and cover overlaps, clipping, strand, unavailable data and source units. Node is needed only for these JavaScript checks, not to run the website.

Availability tests cover below/exactly/above 50%, both strands, precomputed body scores, sparse peaks versus genomic span, missing bounds and the muskrat SIRT6 example. The panel test verifies that both real renderers exclude incomplete selected loci from the fit on raw and log scales while retaining audit values.

Mass-adjusted statistics are checked against independent NumPy full-design OLS and SciPy coefficient t-tests, on every combination of raw/log activity and raw/log MLS. Panel tests check independent toggles, invalid mass exclusions and retention of the gene-body length covariate.
