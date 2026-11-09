# NSC Genome Browser

A Python-backed comparative genome browser for local ATAC-seq data, regulatory annotations, chromatin interactions, and liver RNA/protein measurements. The interface uses plain HTML, CSS, JavaScript, and SVG: no R runtime, Node build step, or cloud service is required.

The current local dataset contains 34 species and 89 ATAC samples. The repository includes the browser, analysis scripts, and a cluster workflow for ATAC-defined candidate super-enhancers. Experimental datasets and generated results are supplied separately.

## Features

- Gene search with exact symbols, unambiguous aliases and supported ortholog mappings, annotation availability counts, and synchronized TSS-relative views.
- Drag to pan tracks; drag species to reorder them; add replicate panels.
- Editable genomic span and signal-bin count, with per-panel or shared signal scales.
- Gene models, ATAC peaks, H3K4me3/H3K27ac regions, human ENCODE4 cCREs, and NMR TE-class filtering.
- NMR predicted Micro-C arcs and a separately labeled hypothalamus ABC test track.
- Human hg38 HepG2 Hi-C heatmap from a local `.mcool` file, automatic resolution, stored-weight balancing when available, and a 700 px heatmap depth.
- Liver transcriptome/proteome scatter plots against MLS or log10(MLS), with Pearson r, p-value, R², and a linear fit.
- Gene and regulatory activity versus lifespan, with optional body-mass correction and phylogenetic regression using the supplied TimeTree covariance.
- Human liver loops, human/mouse dbSUPER tracks, and per-species ATAC candidate super-enhancer shared/union regions.
- Compact single-page vector PDF export in 16:9, 4:3 or full-height layout, per-species annotation toggles, and liver-data CSV export.
- Standalone square regulatory-activity plots with black points, red dashed fits, and BH-adjusted P values with and without body-mass correction.

For required inputs, exact file formats, and setup steps, see the [data requirements guide](docs/DATA_REQUIREMENTS.md).

## Installation and launch

Python 3.9+ is required; Python 3.11 is a reasonable new-environment choice. Availability of binary wheels depends on platform and Python version.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r cross_species_visual_tool/web_app/requirements.txt
python cross_species_visual_tool/web_app/server.py
```

Open http://127.0.0.1:8765/. An alternative port can be supplied with `--port 8766`.

**Data are not included in the code repository.** A fresh clone will not launch successfully until the configured species' bigWigs and transcript indexes are installed. Paths currently follow the local project layout rather than a general configuration file. See [data setup and architecture](docs/CODE_SUMMARY.md).

## Data organization

Reference genomes and annotations are in `reference_genomes/`; experimental inputs are grouped by assay in `data/`. See the [folder guide](docs/FOLDER_LAYOUT.md) for paths and the move manifest.

## Code layout

| Component | Responsibility |
|---|---|
| `cross_species_visual_tool/web_app/server.py` | Local HTTP server, species/sample catalog, interval queries, signal binning, gene anchors, API routes |
| `www/browser.js`, `www/index.html`, `www/style.css` | SVG track rendering, controls, drag interactions, Hi-C display, PDF request generation |
| `chipseq.py`, `ccre.py`, `repeats.py` | Regulatory annotation readers |
| `microc.py`, `abc_tracks.py`, `human_hic.py` | Predicted loops, ABC links, and measured Hi-C matrices |
| `omics.py`, `www/omics.js`, `www/statistics.js` | RNA/protein queries, plotting, and correlation calculations |
| `pdf_export.py` | Sanitized SVG-to-vector-PDF conversion |
| `build_*.py`, `mouse_supplement.py` | Local annotation, TE, lifespan, and workbook preprocessing |
| `test_*.py`, `test_statistics.js` | Backend, data-integration, PDF, and statistical checks |
| `cross_species_visual_tool/plot_regulatory_mls.R` | Regulatory activity/MLS plots and validation of genome-wide BH testing families |
| `pipelines/atac_super_enhancers/` | Filtered-BAM ROSE workflow, SLURM launchers, and replicate-supported region summaries |
| `output/`, `outputs/`, `tmm_pca/` | Analysis source scripts; their experimental inputs and generated results are not included |

## Verification

```bash
python -m pip install -r cross_species_visual_tool/web_app/requirements-dev.txt
cd cross_species_visual_tool/web_app
python -m unittest discover -p 'test_*.py'
node test_statistics.js
```

Node is only needed for the standalone JavaScript test. Many Python tests require the original local dataset; this is not a self-contained demo or a data-free CI test suite.

## Interpretation and deployment

Equal TSS-relative positions and matching gene symbols do not establish orthology. The gene-name resolver distinguishes supported mappings from exact-symbol matches and projected annotations. Signal values retain source units; shared axes do not normalize species. Predicted Micro-C and ABC links are labeled separately from measured Hi-C. Hypothalamus ABC data and HepG2 cell-line data are not matched liver ATAC replicates. Lifespan regressions are exploratory; ordinary and PGLS models are available in the browser. Standalone regulatory plots use ordinary regression.

The server binds to localhost and is intended for local research use. Public hosting would require a separate deployment/security design. Uploading the code to GitHub does not host the Python application or its datasets.

## Project documentation

- [Architecture, data setup, and limitations](docs/CODE_SUMMARY.md)
- [GitHub upload guide](docs/GITHUB_UPLOAD.md)
- [Analysis scripts and regulatory plots](docs/ANALYSIS_SCRIPTS.md)
- [ATAC super-enhancer cluster workflow](pipelines/atac_super_enhancers/README.md)
- [Two research strategies and shared HALPER + IPP framework](docs/cross_species_project_strategies.md)
- [HALPER + IPP workflow](docs/HALPER_IPP_regulatory_workflow.md)

## Licensing and attribution

No software license has been selected yet. Choose one before offering the code for reuse. Third-party data, dependencies, and branding retain their own rights and terms. The Rochester logo in `www/ur-logo.svg` is an existing project asset; inclusion does not grant trademark rights or institutional endorsement.
