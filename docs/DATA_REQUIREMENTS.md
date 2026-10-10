# Data requirements and formats

This guide describes the **current Python NSC Genome Browser**, not a generic importer. Paths below are relative to the project root. GitHub contains code and documentation only; it does not include the biological datasets or generated indexes.

## 1. Minimum data needed to run

For **every species configured in `server.py` → `SPECS`**, provide:

1. At least one ATAC signal file in **bigWig** format, with a filename ending in `.bw`.
2. A gene/transcript index: **`<index-name>_models.tsv.gz` plus its `.tbi` index** in `cross_species_visual_tool/index/`.
3. For Human, **`cross_species_visual_tool/index/Human_anchors.json`**.

The current default configuration has nine species. Hiding a species with the website checkbox does **not** remove its startup requirements. To run a smaller collection, edit `SPECS` to contain only installed species and adjust species-specific readers if introducing new species. There is no upload/import configuration screen.

Raw FASTQ, BAM, reference FASTA, and GTF files are **not read for ordinary browser display** once the bigWigs and transcript indexes exist. FASTA/GTF files are needed for rebuilding transcript indexes. Peak calling and bigWig generation are upstream preprocessing steps; this website does not run alignment or MACS.

### Required ATAC paths and names

| Species ID | Signal folder | Recognized sample prefix | Transcript index basename |
|---|---|---|---|
| Human | `data/atac/selected3/` | Any `.bw` filename | `Human_models.tsv.gz` |
| NMR | `data/atac/mole_rat_peaks/` | `NMR` | `NMR_models.tsv.gz` |
| BMR | `data/atac/mole_rat_peaks/` | `BMR` | `BMR_models.tsv.gz` |
| DMR | `data/atac/mole_rat_peaks/` | `DMR` | `DMR_models.tsv.gz` |
| Mouse | `data/atac/mouse_mm10/` | `wtmice` (case-insensitive) | `Mouse_mm10_models.tsv.gz` |
| Rat | `data/atac/rat/` | `Rat` | `Rat_models.tsv.gz` |
| Macaque | `data/atac/macaque/` | `SRR` | `Macaque_models.tsv.gz` |
| Rabbit | `data/atac/Rabbit/` | `CTR` | `Rabbit_models.tsv.gz` |
| ASM | `data/atac/Africa_spiny_mouse/` | `ASM` | `ASM_models.tsv.gz` |

Examples: `NMR3Liver.cpm.bw`, `WtMice11Liver.cpm.bw`, `122_scaled.bw`. Sample IDs are obtained by removing `.bw`, then `.cpm` or `_scaled`. `NMR1Liver` is deliberately excluded. The old mouse assembly in `data/atac/mouse_legacy/` is not used.

Supply precomputed, appropriately normalized signal. The browser averages supplied values within display bins; it does not compute CPM from raw read counts or normalize abundance across species. Existing human tracks use supplied scaled units; other existing tracks are labeled CPM. For a new collection, update labels to match actual units.

### Transcript index format

This is a **custom eight-column, tab-separated format**, not standard BED8. No header:

| Column | Meaning |
|---|---|
| 1 | Chromosome/scaffold name |
| 2 | Transcript start, zero-based |
| 3 | Transcript end, exclusive |
| 4 | Transcript ID |
| 5 | Gene symbol |
| 6 | Strand: `+` or `-` |
| 7 | Comma-separated absolute exon starts, zero-based |
| 8 | Comma-separated absolute exon ends, exclusive |

Synthetic example (tabs between fields):

```text
chr1	1000	2000	tx_example	GENE_EXAMPLE	+	1000,1800,	1200,2000,
```

Rows must be sorted by chromosome and start, compressed with **BGZF**, and indexed with tabix using BED coordinates. Ordinary gzip alone is not enough for indexed queries. For example, with pysam:

```python
import pysam
pysam.tabix_compress('Human_models.tsv', 'Human_models.tsv.gz', force=True)
pysam.tabix_index('Human_models.tsv.gz', preset='bed', force=True)
```

Human anchors are a JSON array with `symbol`, `chrom`, `strand`, `start`, `end`, `tss`, `transcript`, and `tssMethod`. **Anchor coordinates are one-based inclusive.** Generate this file with `build_human_anchors.py`; it prioritizes MANE Select and other GENCODE transcript tags.

## 2. Reference genomes and rebuilding indexes

Keep each species' FASTA, `.fai`, GTF and any GTF index together in `reference_genomes/`. Assemblies must match the corresponding ATAC data.

The supplied builder handles these five reference sets:

| Reference subfolder | FASTA filename | GTF filename |
|---|---|---|
| `selected3/` | `hg38.fa` | `gencode.v50.primary_assembly.annotation.gtf.gz` |
| `mouse_mm10/` | `mm10.fa` | `Mus_musculus_mm10.ucsc_ncbiRefSeq.gtf.gz` |
| `macaque/` | `Macaca_mulatta.fa` | `Macaca_mulatta.gtf.gz` |
| `Rabbit/` | `Sylvilagus_floridanus.fa` | `Sylvilagus_floridanus.gtf.gz` |
| `Africa_spiny_mouse/` | `Acomys_cahirinus.fa` | `Acomys_cahirinus.gz` (gzipped GTF despite its name) |

Each FASTA needs its matching `<filename>.fai`. The builder reads gzipped GTF text directly; a GTF `.tbi` is not required by this builder. It reads transcript/exon features and transcript identifiers, gene symbols, and strand. It checks bigWig chromosome lengths against the FASTA index.

From the project root:

```bash
mkdir -p cross_species_visual_tool/index cross_species_visual_tool/web_app/data
cd cross_species_visual_tool/web_app
python3 build_species_indexes.py
python3 build_human_anchors.py
```

**Limitation:** the Python builder currently covers Human, Mouse_mm10, Macaque, Rabbit, and ASM only. Supply the existing NMR/BMR/DMR/Rat transcript indexes, or extend the builder for those assemblies. A fresh code clone plus arbitrary FASTA files is not a complete runnable dataset.

## 3. Optional tracks

These inputs add features but are not required for the basic ATAC browser to start.

### ATAC peaks

Place `<sample>_peaks.narrowPeak` or `<sample>_peaks.narrowPeak.gz` beside the matching bigWig. Example: `NMR3Liver.cpm.bw` pairs with `NMR3Liver_peaks.narrowPeak`.

Use standard narrowPeak data with no header. The browser reads the first three tab-separated columns: chromosome, zero-based start, exclusive end. It counts overlapping intervals; it does not call new peaks. ASM peaks are currently disabled because the supplied peak coordinates do not match the signal assembly.

### ChIP-seq regions

Plain BED, at least three columns (`chrom`, `start`, `end`); optional fourth-column name. Use zero-based, half-open coordinates.

- NMR: `data/chipseq/mhetglav3/hetGla-H3K4me3_replicated-peaks_mHetGlaV3.bed` and `hetGla-H3K27Ac_replicated-peaks_mHetGlaV3.bed`.
- Other supported species: `data/chipseq/ChIP_seq/liftover/final/<Species>_<Mark>_<Assembly>.bed`.
- Species/assembly pairs: Human/hg38, Mouse/mm10, Rat/GRCr8, Macaque/T2T-MMU8v2.0. Marks are `H3K4me3` and `H3K27Ac` (case matters).

Use already converted, assembly-compatible intervals. Renaming chromosomes does not perform liftover.

### Human Hi-C matrix

Active file:

```text
data/interactions/human_hic/GSE278978_HepG2-control_merge.mcool
```

Use an HDF5 multi-resolution Cooler file readable by `hictkpy.MultiResFile`, containing resolution groups and a contact matrix at each resolution. The current file has 15 resolutions from 500 bp to 8,192,000 bp; its chromosome lengths match hg38. Its assembly attribute itself is `unknown`, so the compatibility assessment relies on chromosome lengths and the supplied dataset context.

The reader discovers resolutions automatically and selects one targeting approximately 160 bins across the window. A sufficiently coarse resolution is needed to keep chromosome-wide views small. Stored normalization named `weight` is used when present; otherwise raw counts are used. Unsupported chromosomes are reported as unavailable.

To use another matrix, update `human_hic.py` → `PATH`, dataset labels/accession, and the source link in `www/index.html`. The current reader expects `.mcool`; the retained old `.hic` is not active, and changing only the path back to a `.hic` is insufficient.

### NMR predicted Micro-C loops

Folder: `data/interactions/nmr_microc/`.

Names: `NMR3Liver.microc_ch0.enrich1.0.loops.bedpe`, and equivalent NMR4Liver/NMR8Liver files. These follow the selected ATAC sample.

The reader expects **14 tab-separated columns**, in this order:

```text
chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 log2_OE log2_enrichment n_windows_called n_windows_covering
```

Use tabs, not spaces. Both anchors use zero-based, half-open coordinates on mHetGlaV3. A header must begin with `#`; score and enrichment fields must be finite numbers; window counts must be integers. A generic six-column BEDPE alone is insufficient for this reader.

### NMR ABC test links

File: `data/interactions/nmr_abc/mHetGlaV3.primary_hypothalamus_ABC_enhancer.bedpe (1).gz`.

Gzipped, tab-separated text with one header row, then at least eight columns:

```text
chrom1 start1 end1 chrom2 start2 end2 name score
```

The first anchor is the enhancer; the second is the supplied target anchor. The existing target columns represent start-codon coordinates, not necessarily TSS. Target intervals with `start=end` are displayed as one-base points. Gene labels come from the name before the first underscore. This is labeled as a hypothalamus prediction/test track, independent of liver sample selection.

### Human cCREs

File: `data/regulatory_annotations/ENCODE4/encodeCcreRegistry.hg38.bb`.

Use the original ENCODE4 hg38 registry **bigBed schema**, not arbitrary BED converted to bigBed. The reader uses BED column 4 for the name, column 9 for RGB color, and column 10 for category. It is an all-tissue registry, not a liver-specific selection.

### NMR TE repeats

Source: `data/regulatory_annotations/repeats/mHetGlaV3.primary.filteredRepeats.bed.gz`.

Gzipped BED6: `chrom`, `start`, `end`, `class/family`, `score`, `strand`. Class is the part of column 4 before `/`; accepted classes are LINE, SINE, LTR, DNA, Retroposon, and RC. Other classes are filtered out. Run `python3 build_repeats.py` from the application folder. This requires the `NMR3Liver.cpm.bw` reference track and an existing application `data/` directory.

Generated runtime files: `cross_species_visual_tool/web_app/data/nmr_TE.bed.gz` and `.tbi`.

## 4. Optional liver RNA/protein and lifespan data

The plotting API reads these generated files:

- `cross_species_visual_tool/web_app/data/liver_omics.sqlite`
- `cross_species_visual_tool/web_app/data/species_mls.json`

To rebuild them, supply these source workbooks:

| Path | Required identifier columns | Measurement columns |
|---|---|---|
| `data/omics/rna/CrossSpecies_Count.xlsx` | `symbol`, `species`, `sample`, `tissue` | `tpm`, `log.tpm` |
| `data/omics/rna/rnaLongTPM.xlsx` | `symbol`, `species`, `sample`, `tissue` | `tpm`, `log.tpm` |
| `data/omics/protein/protLongDF.xlsx` | `symbol`, `species`, `sample.y`, `tissue` | `quant`, `log.quant`, `quant.norm` |

Column names are case-sensitive. Optional metadata columns are `common`, `source`, and `norm.on`. Rows represent gene/sample measurements. Missing numeric fields become missing values; provide the measurement columns you want to plot. Log values are read as supplied, so their transformation must be known and consistent across sources.

The importer is tailored to the supplied XLSX structure: headers are in the first row of `xl/worksheets/sheet1.xml`, with strings stored in `xl/sharedStrings.xml`. It is **not a general CSV/XLSX importer**; inline-string-only or differently structured workbooks require adapting `build_omics.py`.

Only liver rows are imported. `BMR_Turk_Lu` is excluded as lung. CrossSpecies_Count is primary; rnaLongTPM supplements missing mouse liver records. Scientific species names should be consistent across workbooks and lifespan metadata.

```bash
# From cross_species_visual_tool/web_app
python3 build_omics.py
python3 build_mls.py
```

`build_mls.py` additionally requires:

- `data/omics/lifespan/Anage.zip`, containing `anage_data.txt` with `Genus`, `Species`, and `Maximum longevity (yrs)`.
- A **sibling folder outside this repository**, `../Open4Gene-main/species_MLS.csv` relative to the project root. Required fields are `species`, `status`, `MLS_corrected`, and `source_note`. Adapt this path if necessary.

Lifespans are in years. Verified local reference values take priority over exact scientific-name matches in AnAge. Missing lifespan values are excluded from correlations. Prebuilt `species_mls.json` maps scientific names to `years`, `log10_years`, `source`, and `note` (unavailable species may have `years: null`).

Keep original workbooks available when using the prebuilt SQLite database: it checks workbook sizes and modification times and reports stale sources. If relocating/copying changes timestamps, rebuild the index to refresh provenance.

## 5. Coordinate and compatibility checklist

- All tracks for one species must use the **same assembly**. Human hg19 and hg38, or different mouse assemblies, cannot be overlaid directly.
- Match chromosome/scaffold identifiers and lengths. Some signal/peak lookups handle `chr` prefixes, but annotation readers do not all do so; exact matching is safest.
- BED, narrowPeak, loop anchors, and transcript-index intervals are **zero-based, half-open**. GTF and human-anchor JSON coordinates are **one-based, inclusive**. The browser converts to one-based display coordinates internally.
- Gene symbols anchor comparisons across species; equal symbols/TSS-relative positions do not establish orthologous regulatory elements.
- Do not mix predicted interaction data with measured Hi-C without appropriate labels. The current human matrix is HepG2; the ABC test data are hypothalamus.

## 6. Install, launch, and verify

From the project root, using Python 3.9 or newer:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r cross_species_visual_tool/web_app/requirements-dev.txt
mkdir -p cross_species_visual_tool/index cross_species_visual_tool/web_app/data
# Install the required datasets/indexes described above before launching.
python cross_species_visual_tool/web_app/server.py
```

`requirements-dev.txt` includes runtime dependencies, the XLSX importer dependency, and test dependencies. No R runtime or Node build step is needed. Node is optional for the separate JavaScript statistics test.

Open http://127.0.0.1:8765/. In a second terminal, with the same environment active:

```bash
python -m unittest discover -s cross_species_visual_tool/web_app -p 'test_*.py'
```

The full test suite assumes the original local datasets; it is not suitable unchanged for every reduced/custom collection. Check the species list, search a known gene, verify coordinate ranges, and inspect the optional-track labels. A missing optional file should be treated as unavailable data, not biological absence. Restart the server after replacing files or changing paths to clear cached data.

See [folder layout](FOLDER_LAYOUT.md) and [code summary](CODE_SUMMARY.md) for additional details. Local analysis outputs, old source archives, and the previous `.hic` file are not required to run the active browser.
