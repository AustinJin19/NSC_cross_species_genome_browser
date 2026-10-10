# NSC — Python cross-species ATAC browser

A local Python backend and HTML/CSS/JavaScript interface for the nine species in
this dataset. No R runtime, Node build, external fonts, or cloud service is used.

## Run

Requires Python 3.9+ and the packages listed in `requirements.txt`:

```sh
python3 -m pip install -r requirements.txt
python3 server.py
```

Open http://127.0.0.1:8765. You can launch `server.py` by absolute path from any
working directory. Use `--port 8766` to choose another port. It binds only to
127.0.0.1 and is intended for local research use, not a public deployment.

## Use

- Search a gene symbol, starting with HAS2. Autocomplete searches the local indexes.
- Show or hide any species and select its individual liver replicate.
- Pan upstream/downstream together, zoom, or enter any positive integer span and choose 10–5,000 bins.
- Hover over signal to see a genomic bin center and its mean CPM value.
- Switch between independent and shared signal ranges.
- Add another panel for any species to compare replicates simultaneously; remove extra panels independently. Hidden species retain their panel settings.
- Set minimum and maximum CPM per panel, or leave either blank for automatic limits. Set selects independent scaling; shared mode temporarily uses a common automatic range. Values outside manual limits are visually clipped; exports retain original values.
- Gene search reports annotated samples and species across the full catalog, including in autocomplete. Counts refer to annotation availability, not measured accessibility; all samples of a species share its annotation.
- Toggle narrowPeak calls and gene models. The legacy signal CSV API remains, but the signal-export button has been removed.
- Export PDF downloads a vector rendering of the current panels, scales, sample selections, peaks, and visible gene annotations. Multi-page exports repeat the header and ruler, with four displayed panels per page.

## Data and coordinate conventions

Paths are resolved relative to this project, without the stale absolute paths in
its original R config. The backend reads existing `index/*_models.tsv.gz` files
and their tabix indexes with pysam, and local bigWigs with pyBigWig. It does not
read RDS files or invoke R. NMR1Liver is excluded from the browser catalog; NMR3Liver is the default naked mole-rat sample. Source files are retained.
Replicates are discovered in the mouse, rat, and
mole_rat_peaks folders; each matching narrowPeak file is loaded on first use.
The original data and original R application are unchanged.

Transcript indexes have nine tab-delimited columns: chromosome, 0-based start,
exclusive end, transcript ID, gene symbol, strand, comma-separated exon starts,
comma-separated exon ends, and biotype. Exon starts are 0-based. Displayed
coordinates are converted to 1-based inclusive coordinates.

Symbols are matched case-insensitively. For nonhuman species, the merged gene boundary is
used as the TSS (highest end on the minus strand). Human uses the prioritized transcript described below. Signal and coordinates are
reversed on minus-strand genes so upstream appears to the left in all panels.
The requested number of mean bins is returned (700 by default), capped at the number of bases; uncovered positions inside a bigWig sequence
are treated as zero, whereas positions outside the sequence are missing.
A window labeled 10 kb spans -5 kb to +5 kb, including the center base.

The CSV includes panel ID, sample, species, gene, strand, relative position from TSS,
genomic bin center, and mean CPM; bin centers can be fractional. The shared
scale uses the same numeric range, without performing additional normalization.

## Current limitations

This first version uses symbol matching, not an orthology database, sequence
alignment, or liftOver. For example, TP53 and TRP53 are different search keys.
Equal distance from a TSS does not imply homologous regulatory sequence.
Gene records sharing a symbol, chromosome, and strand are collapsed together;
when multiple chromosome/strand matches remain, the longest span is chosen.
Assembly accession metadata is not available in the supplied index and is not
inferred from species names.

The longest transcript represents each gene. At wide windows at most five gene
models are displayed, prioritizing the searched gene and reporting the total.
There is no arbitrary coordinate entry, reference-base track, alternate-TSS
selector, or multi-replicate aggregation yet. BigWig values are queried using exact per-bin sums; no fixed upper window-span limit is imposed.

## Verification

```sh
python3 -m unittest discover -s . -p 'test_*.py' -v
```

The integration tests use the local dataset to check both strand orientations
against raw bigWig bases, synchronized panning, overlapping peaks, missing genes,
out-of-sequence behavior, replicate changes, and wide-window gene selection.

The sidebar reuses `ur-logo.svg` from `Open4Gene-main/website_NSC`, paired with the Upstate Nathan Shock Center name. PDF export uses ReportLab and svglib; install requirements before running. The server retains only the eight most recent PDF downloads in memory.

### Liver multi-omics and lifespan correlations

The liver section follows the main gene search. RNA comes from
`data/omics/rna/CrossSpecies_Count.xlsx`; protein comes from `data/omics/protein/protLongDF.xlsx` under `data/omics/` in the project
root. Rebuild the read-only SQLite index with `python3 build_omics.py` (requires
lxml). Workbook dimensions are unreliable, so the importer streams actual XML
rows. Workbooks are never modified. Only rows labeled Liver are imported;
BMR_Turk_Lu is explicitly excluded as Lung, while BMR_Turk_Liv is the liver sample.
Mouse liver records missing from the primary source are supplemented from
`data/omics/rna/rnaLongTPM.xlsx`, with per-record workbook provenance. Comparing 928,299 shared
liver records showed exactly identical TPM and log.tpm. Existing primary records
take precedence; the supplement is restricted to Mus musculus liver.

RNA has independent measurement controls: supplied log.tpm (default) or TPM.
The RNA workbook has no quant.norm column, so normalized RNA values are missing,
not filled from the old workbook’s quant.norm column. Protein retains normalized, log and raw choices.
Each point is a species median. OLS with an intercept, Pearson r, two-sided p-value and R² use species
as observations, requiring at least three species and variation on both axes.
Choose MLS in years or log10(MLS). Fits are exploratory, without phylogeny or study
adjustment; absent measurements remain absent.

`data/species_mls.json` records local lifespan values and provenance. Rebuild with
`python3 build_mls.py`, using the sibling Open4Gene-main/species_MLS.csv verified
MLS_corrected entries, then exact-name matches from project data/omics/lifespan/Anage.zip. Unmatched
species are excluded from correlations. These are source snapshots, not current
longevity-record claims.

Export liver data downloads all liver records for the gene, including workbook,
source row, measurement values, lifespan and its provenance, regardless of the
plot species filter. The main Export PDF button exports ATAC panels only.

Validation: `python3 -m unittest discover -p 'test_*.py'` and
`node test_statistics.js` from this directory.

### NMR ChIP-seq regions

The NMR panels include optional H3K4me3 and H3K27ac lanes from the replicated-peak
BED files in project `data/chipseq/mhetglav3/`. These are region annotations, not quantitative
ChIP signal. Coordinates are labeled mHetGlaV3 in the source filenames; all input
intervals fit the chromosome names and sizes in the NMR bigWig. BED starts are
converted from 0-based to 1-based for display. The lanes follow the same gene
orientation, pan, and zoom as ATAC, appear in each NMR panel, and are included in
PDF exports when enabled. Sidebar checkboxes independently control the two marks.
No tissue or replicate match to the ATAC samples is inferred from these files.

### NMR TE repeats

The NMR-only TE lane uses `data/regulatory_annotations/repeats/mHetGlaV3.primary.filteredRepeats.bed.gz` from the project
root. Run `python3 build_repeats.py` to rebuild the sorted BGZF/tabix index in data/.
The source BED is unchanged. Classes LINE, SINE, LTR, DNA, Retroposon, and RC are
included; simple repeats, low complexity, satellites, RNA annotations, Unknown,
and artefacts are not classified as TEs for this lane. Unsupported fragment
sequences and invalid/zero-length intervals are excluded. The build metadata is
`data/nmr_TE.bed.json`: 3,575,155 indexed regions, 47,533 TE rows on unsupported
sequences and 1,418 invalid TE intervals excluded. Coordinates on supported
sequences were checked against the NMR bigWig chromosome bounds.

The display converts BED starts to 1-based inclusive coordinates, follows the
panel's strand and window, and retains the source class/family, strand and score
in hover details. Each TE class occupies its own colored subrow, in sidebar
legend order. The sidebar toggle controls all NMR panels, including PDF export.

### Additional ATAC species and mm10 mouse (2026-09-30)

Before adding human, the browser loaded eight species and 23 samples. Mouse uses only the three
`data/atac/mouse_mm10/` bigWigs, their matching peaks, and the mm10 RefSeq annotation.
The old mouse files/index are preserved but are no longer used by the web app.
Rhesus macaque (`Macaca mulatta`, four SRR tracks), eastern cottontail
(`Sylvilagus floridanus`, two CTR tracks), and African spiny mouse
(`Acomys cahirinus`, two ASM tracks) are loaded from the project-root folders.

Rebuild the four transcript indexes with
`python3 cross_species_visual_tool/web_app/build_species_indexes.py`.
The builder checks bigWig chromosome lengths against the supplied FASTA indexes,
extracts bed2gtf gene symbols from transcript IDs, and maps spiny-mouse scaffold
names using the explicit accession/scaffold pairs in the supplied FASTA headers.
Index provenance is recorded in adjacent `.json` files. The input files are unchanged.

Spiny-mouse peak calls are disabled: the supplied narrowPeak files use CM057
chromosomes, whereas its bigWigs and mapped annotation use JAULSH scaffolds.
A compatible peak file or validated coordinate conversion is required to enable them.
RNA/protein measurements retain their existing workbook sources; adding ATAC
species does not create missing RNA/protein observations.

### Human hg38 panel (selected3)

The human panel uses `data/atac/selected3/122_scaled.bw`, `151_scaled.bw`, and
`152_scaled.bw`, their matching gzip-compressed narrowPeak files, and GENCODE
v50 primary-assembly GTF. The supplied hg38 FASTA index validates bigWig sequence
lengths. Transcript features on sequences absent from the signal tracks are
omitted during indexing. Sample labels are 122, 151, and 152. Human signal is
labeled "scaled", preserving the supplied values without assuming CPM units.
The browser contains nine species and 26 samples. Signal CSV responses use
`mean_signal` plus `signal_unit` to distinguish the supplied scales.

Human TSS anchors use transcript tags from the supplied GENCODE GTF, in priority
order: MANE Select, Ensembl canonical, APPRIS principal 1, basic, then longest
transcript. Ties use transcript length and ID. The selected transcript is also
preferred in the gene-model lane, and its ID/selection method appear in the panel.
This avoids anchoring to the outermost boundary of all alternative isoforms.
Rebuild anchors with `python3 cross_species_visual_tool/web_app/build_human_anchors.py`
(or rebuild the human index). ALDH1A2 uses MANE ENST00000249750.9 at
chr15:58,065,711 on the minus strand (hg38, 1-based).

### NMR predicted Micro-C loops

The NMR panel includes an independently toggleable arc lane from project-root
`data/interactions/nmr_microc/NMR{3,4,8}Liver.microc_ch0.enrich1.0.loops.bedpe`, matched to the selected ATAC sample.
These are labeled EPCOTv2 predictions on mHetGlaV3, based on the user's confirmation.
All 156,288 loop records pass coordinate-bound checks against the NMR signal assembly.
BEDPE input is 0-based half-open; the API and tooltips use 1-based inclusive anchors.
Loops are included when either anchor overlaps the view. Off-screen partners use
 dashed arcs ending at the view boundary. Arc height encodes displayed separation,
not interaction strength. Hover reports both original anchors, log2 O/E, log2
 enrichment, and called/covering window counts. At most 300 loops per view are
shown, ranked by score, with the full matching count disclosed. PDF export includes
visible arcs. No source loop files are modified.

### Hypothalamus ABC test track

The NMR ABC test lane reads the project-root gzipped hypothalamus ABC BEDPE
(29,185 links), independently of the selected liver sample. Purple arcs link the
source enhancer interval to the supplied gene anchor. The header calls this
anchor `startCodonStart/End`, not TSS. Equal start/end values are interpreted as
zero-based points and displayed as one-base intervals; original anchor values
are retained in hover details. This coordinate convention should be confirmed
against the generating pipeline before quantitative analysis. All anchors fit
the NMR chromosome bounds. ABC score, gene ID, and annotation are available on
hover. Either-anchor overlap, a disclosed top-300 score cap, dashed off-screen
partners, strand orientation, and PDF inclusion follow the Micro-C lane.

The active Micro-C dataset uses enrich1.0 (`min_enrich=1.0` in the supplied metadata). Original files are retained. Sample-specific tracks remain matched to NMR3/4/8; consensus files are not substituted for individual samples.

### Human Hi-C heatmap

`human_hic.py` reads `data/interactions/human_hic/GSE278978_HepG2-control_merge.mcool` using hictkpy. This is
The supplied HepG2 control merged matrix has chromosome lengths matching hg38; it is not a matched primary-liver
ATAC replicate. A 700 px triangular heatmap follows the human genomic window and
strand. Resolution is selected automatically to bound matrix size. weight-balanced
contacts are shown when available; fallback raw counts are explicitly labeled.
Color uses log1p scaling with a view-specific 98th-percentile cap. Masked values
are gray. The heatmap can be toggled and is included in the track PDF.

See the repository-level README and `docs/CODE_SUMMARY.md` for the current
architecture and source-only GitHub distribution requirements.
