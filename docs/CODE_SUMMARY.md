# Code summary and data setup

## Architecture

The Python `ThreadingHTTPServer` serves the static interface and JSON APIs. At startup, `load_catalog()` discovers configured bigWigs and loads transcript indexes. A browser query resolves the gene in each species, chooses its anchor, computes a strand-aware genomic interval, and reads exact binned bigWig sums. The frontend renders the returned values and annotation intervals as SVG. Both axes of the human Hi-C triangle follow the human panel's genomic orientation.

Human anchors prioritize MANE Select and other transcript tags. Other species use the merged gene boundaries. Ambiguous gene symbols select the longest locus. BED/tabix inputs are generally zero-based, half-open; displayed intervals are one-based inclusive. Outside-chromosome signal is missing, while uncovered positions inside the signal chromosome are treated as zero.

The Hi-C reader queries the local matrix at a resolution chosen to bound allocation to approximately 160 bins per axis. Color uses log1p scaling, capped at the view's 98th percentile. Missing/masked values remain distinct from zeros. Each request uses its own reader; a small interval cache avoids repeat work.

## API routes

| Route | Purpose |
|---|---|
| `GET /api/catalog` | Species, sample labels, and defaults |
| `GET /api/genes?q=...` | Symbol completion and annotation availability |
| `POST /api/browse` | Gene, half-window, pan offset, bins, and requested panels → track data |
| `GET /api/omics?gene=...` | Liver RNA/protein records and metadata |
| `GET /api/omics-export?gene=...` | Liver records as CSV |
| `POST /api/export-pdf` | Current SVG panels → temporary PDF download URL |
| `GET /api/download-pdf/<token>` | Download a generated PDF |
| `GET /api/export?view=...` | Legacy signal CSV endpoint; no visible export-signal button |

## Local data layout

Preserve these paths relative to the repository root, or deliberately update the corresponding reader configuration. The source-only archive excludes all of these datasets and generated indexes.

| Location | Contents / consumer |
|---|---|
| `cross_species_visual_tool/index/` | `*_models.tsv.gz` plus `.tbi` indexes; `Human_anchors.json` |
| `cross_species_visual_tool/mole_rat_peaks/` | NMR, BMR, DMR bigWigs and narrowPeak files |
| `cross_species_visual_tool/rat/` | Rat signal and peaks |
| `selected3/` | Human scaled bigWigs, peaks, hg38 annotations |
| `mouse_mm10/`, `macaque/`, `Rabbit/`, `Africa_spiny_mouse/` | Additional species' signal, peaks, and annotation inputs |
| `mhetglav3/` | NMR histone-region BED files |
| `ChIP_seq/liftover/final/` | Other species' assembly-compatible histone peaks |
| `ENCODE4/encodeCcreRegistry.hg38.bb` | Human cCRE registry |
| `loops/` | Sample-specific NMR enrich1.0 BEDPE predictions |
| `4DNFICSTCJQZ.hic` | Human HepG2 hg38 Hi-C; [4DN source](https://data.4dnucleome.org/files-processed/4DNFICSTCJQZ/) |
| `mHetGlaV3.primary_hypothalamus_ABC_enhancer.bedpe (1).gz` | NMR ABC test input |
| `mHetGlaV3.primary.filteredRepeats.bed.gz` | Source repeat annotations |
| `CrossSpecies_Count.xlsx`, `rnaLongTPM.xlsx`, `protLongDF.xlsx` | Primary RNA, mouse RNA supplement, and protein workbooks |
| `cross_species_visual_tool/web_app/data/` | Generated omics SQLite, lifespan JSON, and TE BGZF/tabix files |

Startup currently requires signal and transcript indexes for all configured species, even when some are hidden in the UI. There is no automatic downloader or bundled synthetic dataset. Configured file names can be inspected in `server.py`, `chipseq.py`, and the respective readers.

## Preparing indexes

Run builders from `cross_species_visual_tool/web_app` after restoring source inputs:

```bash
python build_species_indexes.py
python build_human_anchors.py
python build_repeats.py
python build_omics.py
python build_mls.py
```

These commands are not a universal end-to-end importer. `build_species_indexes.py` covers human, mm10 mouse, macaque, rabbit, and ASM; the original NMR/BMR/DMR/rat transcript indexes must be supplied separately. Its exact GTF/FASTA-index filenames are defined in `SOURCES`. The omics builder requires lxml and the workbook sources; the lifespan builder requires the omics database, project `Anage.zip`, and sibling `Open4Gene-main/species_MLS.csv`.

## Important dataset decisions

- NMR1Liver is excluded. NMR3/4/8 predictions follow the selected sample.
- Human bigWigs retain supplied scaled units; mouse uses mm10 rather than the old assembly.
- ASM peak calls are disabled because their chromosome names do not match the signal assembly.
- BMR_Turk_Lu is excluded from liver RNA because it is lung. Primary RNA records take precedence; missing mouse liver records are supplemented from rnaLongTPM.
- Loop arcs show a disclosed top-300 subset when necessary. Hi-C is a matrix, not a list of statistically called loops.
- PDF exports contain four displayed panels per page, including duplicated species panels. Liver scatter plots are not part of the track PDF.

## Scope

The current implementation visualizes supplied data and predictions. It does not train EPCOTv2, run HALPER/IPP, build a Cactus alignment, or implement the proposed multispecies accessibility model. Those plans are documented separately.
