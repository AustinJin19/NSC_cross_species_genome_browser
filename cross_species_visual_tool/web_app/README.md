# NSC — Python cross-species ATAC browser

## Gene names and ortholog searches (30 October 2026)

The search, gene activity, regulatory activity and liver omics endpoints share
`gene_names.py`. The portable reference is `data/gene_names.json.gz`; deploy it
with the code. Original annotation and source symbols are retained. Searching
either `HBE1` or `Hbb-y` resolves mouse `Hbb-y`. Rabbit and African spiny mouse
use their explicitly mouse-derived `ENSMUST...#Hbb-y#...` annotations, labelled
as projected annotations rather than independently verified orthologs.

HGNC aliases are used only when they identify one approved symbol. Approved
human symbols take precedence over synonyms. HCOP mappings cover mouse, rat
and macaque and require at least three supporting resources, including
Ensembl, NCBI or HGNC, plus a unique reciprocal partner within that filtered
mapping. Ambiguous mappings are flagged and not pooled. An exact local symbol
has priority; matching symbols in other species do not by themselves establish
orthology. Multiple loci of the same resolved symbol retain the existing
longest-locus browser / highest-median activity selection rules.

RNA/protein matching is performed independently per assay and does not infer
mouse-derived annotation provenance from the genome annotation. CSV exports
include the original gene symbol and matching method. Existing analysis
exports and precomputed databases are not rewritten.

Sources: [HGNC](https://www.genenames.org/download/statistics-and-files/) and
[HCOP](https://www.genenames.org/help/hcop/). The bundled snapshot records source
file checksums and download URLs. Rebuild with `build_gene_names.py --sources
SOURCE_DIRECTORY --hgnc HGNC_COMPLETE_SET.tsv`; the directory must contain the
three human_mouse, human_rat and human_macaque HCOP fifteen-column `.txt.gz`
files. Tests: `python3 -m unittest test_gene_names` from this directory.

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
- Use **Hide annotations / Show annotations** in each species panel to control its gene track independently. The sidebar **Gene annotations** checkbox is the global switch; turning it back on restores each panel's setting. Panel settings persist while changing genes, samples, or the displayed window, and PDF exports reflect the visible tracks.
- PDF exports trim empty bands between tracks and unused space below annotations, with compact species separators. Complete shapes move vertically without changing their genomic x positions, signal amplitudes, or loop geometry; the interactive browser layout stays unchanged.
- Export PDF downloads all selected panels as vectors on one page. The **PDF page** selector offers **Slide 16:9** (960 × 540 pt), **Slide 4:3** (720 × 540 pt), and **Full height**. Both slide formats use one vertical stack in browser order, with labels on the left, full-width plots on the right, and one shared ruler. Track geometry is compressed vertically to fit the page; text glyphs retain their proportions and genomic positions. Slide exports have a white background and no title block or footer. Full height grows vertically and preserves plot aspect ratios for detailed views. Dense gene labels are smaller on fixed slides; no tracks or annotations are removed. Human liver loop metadata is stripped during SVG validation while the loop arcs are retained.

## Data and coordinate conventions

Paths are resolved relative to this project, without the stale absolute paths in
its original R config. The backend reads existing `index/*_models.tsv.gz` files
and their tabix indexes with pysam, and local bigWigs with pyBigWig. It does not
read RDS files or invoke R. The website uses the pre-QC sample set: 89 samples
across 34 species, including the low-TSS samples and all three bowhead samples.
The TSS enrichment filter in `data/sample_qc_exclusions.json` is disabled
(`enabled: false`) at the user's request. Its candidate list is retained for
reference. `always_exclude` preserves the earlier NMR1Liver exclusion; that sample
was already absent before the TSS filter was introduced. NMR3Liver remains the
default naked mole-rat sample. Gene activity, peak counts and CPM-in-peaks all
use the restored sample set and recompute species medians before choosing the
highest-activity locus. QC-filtered analyses remain separate historical outputs;
they are not the website's current results. Original measurements are unchanged.
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

Spiny-mouse peak calls use corrected narrowPeak files in
`data/atac/liver_atac_bundle/Acomys_cahirinus/`, matching the JAULSH scaffolds in
its bigWigs and annotation. The old incompatible CM057 peak files were removed.
RNA/protein measurements retain their existing workbook sources; adding ATAC
species does not create missing RNA/protein observations.

### Human hg38 panel (SRR CPM tracks)

The human panel uses four paired bigWig/narrowPeak samples in
`data/atac/Homo_sapiens/`: SRR28763537, SRR28763586, SRR28763842, and
SRR28763939. SRR28763586 is the default. The `.cpm.bw` files retain their supplied
CPM values; no additional normalization is applied. The import manifest in that
folder records original paths, SHA-256 checksums, sizes and peak counts. All peak
intervals are within their paired bigWig chromosome bounds. All 25 exact-name
matches to the supplied hg38 FASTA index agree in length, and every contig in the
existing human annotation index is present. Extra UCSC contig aliases have not
been independently sequence-verified. Raw-read processing metadata was not supplied.

GENCODE v50 and the existing hg38 gene models/TSS anchors are reused; the reference
files remain in `reference_genomes/selected3/`. The old scaled tracks are retained
in `data/atac/selected3/` for historical analyses and are not active samples.
Live gene-activity and regulatory queries reconcile the sample set with the current
catalog and recalculate scores on demand when it differs from the stored matrix.
Human integrated activity is now CPM·bp. Human is excluded from the live lifespan
plots and fits, with values retained in the tables. Macaque is included when eligible,
with the same region, measurement and covariate requirements as other species. Historical matrices and analysis
exports are not rewritten by this panel update. Signal CSV responses include
`mean_signal` and `signal_unit`.

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

### NMR hypothalamus Hi-C

The NMR Hi-C heatmap replaces the former ABC test track and reads
`data/interactions/nmr_microc/NMR_hypothalamus_merged.mcool` on mHetGlaV3.
This hypothalamus matrix is shared across the NMR replicate panels; it is not
liver-matched data. The browser selects among the supplied 5 kb–1 Mb resolutions
to keep regional queries bounded, uses balancing weights when available, and
preserves masked bins. Toggle it independently using **NMR Hi-C heatmap**.
The existing predicted Micro-C loops and human Hi-C heatmap remain separate.

### Mouse liver ZT10 Hi-C

The mouse mm10 panels read
`data/interactions/mouse_hic/mouse_liver_ZT10_merged.mcool`.
Enable **Mouse Hi-C heatmap** to show the merged liver ZT10 matrix beneath
the other mouse annotations. The same matrix is used for each mouse ATAC
replicate; it is not a replicate-specific contact map. The heatmap uses the
shared reader, stored balancing weights, automatic 5 kb–1 Mb resolution,
masked-bin handling, and 700 px depth. It is included in PDF exports.
The mcool data file must be copied separately when deploying the code.

## Gene and regulatory activity

The browser includes separate gene-activity and regulatory-activity panels with MLS correlations. Regulatory activity offers peak count and CPM signal integrated only within peaks for TSS ±20/50/100/200/500 kb windows. See [GENE_ACTIVITY.md](GENE_ACTIVITY.md) for definitions, covariates, required generated files, and deployment steps. Hi-C heatmaps are off by default and can be enabled in Track display.

### Phylogenetic generalized least squares (PGLS)

Both activity panels have a **Regression** selector: ordinary regression,
PGLS with Pagel's λ estimated by maximum likelihood in [0, 1], or PGLS with
λ fixed at 1 (Brownian motion). The response is the selected raw or log10
species-median activity; the predictor is the selected raw or log10 MLS.
Gene-body models retain log10(annotated gene length in kb), and the existing
body-mass checkbox adds log10(body mass in g). Human is excluded and macaque is eligible; bowhead-option,
MLS, coverage and activity filters are identical to the ordinary model.
Species absent from the tree are explicitly excluded in PGLS only.

The residual covariance matrix uses shared evolutionary time:
`C[i,j] = root-to-MRCA(i,j) branch length`, with root-to-tip length on the
diagonal. It is not a matrix of pairwise tip distances. For λ models,
off-diagonal entries are multiplied by λ while diagonal entries stay fixed.
Each fit uses the principal submatrix for the included species in the exact
data order, preserving the original tree root. The displayed line holds
covariates at their arithmetic means across included species.

The bundled tree and covariance come from the
[TimeTree 6 mammal tree](https://timetree.org/downloads/biodiversity/TT6_Mammals.nwk).
They contain the original 34-species analysis tree plus bowhead whale, added
from the same source, covering all 34 species currently in the browser.
Source SHA-256, branch units, taxon names and build provenance are recorded in
`data/phylogenetic_covariance.json`. The derived Newick is
`data/timetree_pgls.nwk`; `www/phylogeny.js` bundles the covariance for the UI.
Rebuild with `python3 build_phylogeny.py` (Biopython and NumPy required for
this build step only). The build validates the matrix and checks that all
original 34-species branch distances remain unchanged.

`www/pgls.js` fits GLS by Cholesky whitening and QR decomposition. It reports
the MLS coefficient, standard error, two-sided coefficient t-test p-value,
λ, GLS R², sample size and residual degrees of freedom. The t test is
conditional on λ; estimated-λ p-values use the standard plug-in
approximation. GLS R² compares to an intercept-only model at the same λ;
it is not Pearson r². P-values are nominal and do not update earlier
genome-wide FDR lists. Invalid or singular fits do not fall back to OLS.

The browser needs no additional runtime dependency. When deploying, include
`www/pgls.js` and `www/phylogeny.js` and restart `server.py` so their new
static routes are available. Numerical fixtures were independently generated
with R's [caper::pgls](https://search.r-project.org/CRAN/refmans/caper/html/pgls.html).
Run `node test_pgls.js`, `node test_statistics.js` and
`node test_activity_panel.js`. To regenerate the independent fixtures, run
`Rscript generate_pgls_fixtures.R` with ape, caper and jsonlite installed.
# ATAC candidate super-enhancers: shared and union regions

The **Show ATAC super-enhancers** checkbox enables species-level candidate SE regions.
Use **ATAC super-enhancer regions** to select **Shared in ≥2 samples** from
`species_summary/*_candidate_se_shared_ge2.bed`, or **Union (≥1 sample)** from
`species_summary/*_candidate_se_union.bed`. Union includes candidate SE calls
from any sample, including species with only one sample; shared regions require
at least two distinct samples at every base. Adjacent qualifying segments are
merged; supporting sample identities may change along a region. The track
is independent of the selected ATAC replicate and the human/mouse dbSUPER tracks.
These ATAC candidates do not establish tissue specificity or functional validation.

Install or update from the project root:

```bash
python3 cross_species_visual_tool/web_app/import_atac_super_enhancers.py
```

The importer verifies sample identities, exact support-table agreement, counts,
and chromosome names/bounds against all contributing browser bigWigs before
copying the BED files unchanged into `data/regulatory_annotations/atac_super_enhancers/`.
`metadata.json` records sample IDs, counts, checksums, and unavailable species.
Both region sets are checked against the source support table. BED input is
0-based half-open; the browser converts to 1-based inclusive. Deploy this data
directory along with the code. The sets remain separate, with shared selected
by default. Single-sample shared regions and absent files are explicitly identified.
These display tracks do not change activity or lifespan calculations.

# Human liver DeepLoop loops (hg38)

The **Human liver Hi-C loops** checkbox shows GSE167200 liver contact pairs in
the human panel, independently of the HepG2 Hi-C heatmap. The loop track is on
by default. These are the published top300K ranked DeepLoop contact pixels,
not 300,000 independently validated or FDR-filtered loops. Hover over an arc
for the original DeepLoop score, hg38 coordinates and original hg19 anchors.
The same tissue resource is shown for each human ATAC replicate.

Installed data: `data/interactions/human_liver_deeploop/GSE167200_Liver.top300K.hg38.bedpe.gz`
(paths in this section are relative to the project root). BEDPE uses 0-based
starts and exclusive ends; browser coordinates are 1-based inclusive.
Both source anchors are converted with UCSC hg19ToHg38 liftOver at minMatch
0.95. Source bins are interpreted as 1-based inclusive from their contiguous
boundaries and minimum start of 1. Pairs are rejected if either anchor fails
mapping, maps ambiguously, exceeds target chromosome bounds, changes length
by more than twofold, or maps a source cis pair to different chromosomes.
The current conversion retains 299,651 of 300,000 pairs. The accompanying
`conversion_summary.json`, `anchor_mapping.tsv.gz` and `rejected_pairs.tsv.gz`
record provenance, coordinate conventions and exclusions. The original input
is preserved. The browser queries either anchor, deduplicates pairs, and draws
up to 300 highest-scoring overlapping pairs. Dashed arcs indicate off-screen
anchor midpoints. No genomic-distance or gene-specific filtering is applied.

Rebuild from the project root:

```bash
python3 cross_species_visual_tool/web_app/build_human_liver_loops.py
```

Use `--liftover`, `--chain`, `--source`, `--sizes` and `--output` to override the
local defaults, including a platform-appropriate UCSC executable. Only the
converted BEDPE file is needed to serve the track; copy it with the rest of
`data/` when deploying. Source: [GSE167200](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167200).
