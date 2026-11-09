# Cross-species ATAC-seq browser

Whole-genome browser that stacks one panel per species on a shared axis measured
as distance from the TSS, so the same locus can be compared across assemblies.

Search a gene to anchor every panel on that species' own TSS, or type a region
for a single species to move that panel independently.

## Files

| File | Purpose |
|---|---|
| `config.tsv` | Points at your FASTA, GTF and bigWig per species |
| `build_index.R` | One-time preprocessing; reads the GTF once |
| `app.R` | The browser |
| `atac-browser-demo.html` | Standalone prototype with simulated data, no R needed |
| `app_upload_version.R` | Earlier upload-based version, kept for reference |

## Setup

```r
install.packages(c("shiny", "ggplot2", "data.table"))
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges", "Rsamtools", "Biostrings"))
```

ggplot2 3.4+ is required (`linewidth`).

Edit `config.tsv`, then:

```sh
Rscript build_index.R config.tsv index/
```

```r
shiny::runApp("app.R")     # or: ATAC_INDEX_DIR=/path/to/index R -e 'shiny::runApp("app.R")'
```

## Why there is a build step

A GENCODE human GTF is around 1.6M lines. `rtracklayer::import()` on it takes
1–3 minutes and several GB of RAM. Doing that per session makes the app
unusable, so `build_index.R` reads each GTF exactly once and writes:

- `<species>_genes.rds` — one row per gene, held in memory for search and TSS
  lookup. About 3 MB for human.
- `<species>_models.tsv.gz` + `.tbi` — transcript models, bgzipped and tabix
  indexed. Range queries stay off the heap, which matters with several species
  open at once.
- `<species>_meta.rds` — file paths, chromosome sizes, and the sequence-name
  maps described below.

bigWig and FASTA are left alone; both are already indexed formats with fast
random access. The FASTA gets a `.fai` written next to it if one is missing,
so that directory must be writable once.

Rebuild the index when your annotation or bigWig files change. Nothing else
caches.

## Sequence naming

The most common cause of empty panels is a GTF that says `1`, a bigWig that
says `chr1`, and a FASTA header that says `1 dna:chromosome ...`. The build
step resolves this once and stores explicit maps, rather than guessing at
query time. Watch the build output:

```
Human
    reading GTF (this is the slow step)
    reconciling sequence names
    bigWig   all 194 sequences matched
    FASTA    189 of 194 matched; unmatched: KI270302.1, KI270303.1, ...
    dropping 5 sequence(s) absent from bigWig or FASTA
```

Sequences missing from either the bigWig or the FASTA are dropped from the
index, because a gene on them could never be drawn. If a whole assembly shows
zero matches, the naming conventions differ by more than a `chr` prefix and
you will need to rename one of the files with something like
`bigWigToBedGraph` plus a rename step, or `seqkit replace` on the FASTA.

## The coordinate model

```
genomic position = anchor[species] + nudge[species] + pan + relative_x
```

- **anchor** — per species. Set by gene search (that species' canonical TSS,
  the most upstream transcript start) or by typing a region, in which case it
  is the region midpoint.
- **nudge** — per species, in bp. For when the annotated TSS is wrong in one
  assembly and you want to line the promoter peaks up by eye.
- **pan** — shared. Every panel moves together, and it moves in the anchor's
  strand direction so panels stay aligned when a gene is on opposite strands
  in different species.
- **window** — shared half-width.

Typing a region with a span also sets the shared window to match, which is
almost always what is meant by pasting coordinates.

`Recentre` clears pan and all nudges. It does not clear manual anchors —
search a gene again to reset those.

## Scaling

Signal is normalised to 0–1 within each panel, and the true maximum is printed
on the panel. **Independent** divides each panel by its own maximum, so you are
comparing peak shape. **Shared** divides everything by the global maximum, so
you are comparing magnitude — only meaningful once the per-species scale
factors put the tracks on one footing.

Compute those factors upstream. `bamCoverage --normalizeUsing RPGC
--effectiveGenomeSize <n>` is the usual choice for cross-species work, because
it corrects for the mappable genome size differences that make raw CPM
misleading between, say, human and zebrafish. Enter the resulting factor in the
box next to each species.

## Tracks

- **Signal** — binned mean from the bigWig. The bin count is shared; at wide
  windows increase it or narrow peaks will be averaged away.
- **Gene models** — collapsed to the longest transcript per gene by default.
  `All transcripts` shows up to 8, longest first.
- **Sequence / GC** — GC fraction in bins below 500 kb of total view, and
  individual bases below 400 bp. Both are reverse-complemented for
  minus-strand anchors, so the sequence shown reads in the direction of
  transcription.

## Performance

Query cost is dominated by the bigWig read, which is proportional to window
width rather than genome size. A ±5 kb view across 5 species is a few tens of
milliseconds. A ±5 Mb view is seconds, mostly spent decompressing. If wide
views are a common use case, generate reduced-resolution bigWigs and pick
between them by window size.

Memory at rest is one gene table per species, so roughly 3–5 MB each. The
tabix models and the bigWigs are read on demand and never held.

## Known limits

- Cross-species gene matching is by upper-cased symbol. This is fine among
  mammals and degrades outside them. The earlier ortholog-table approach in
  `app_upload_version.R` is the fix; wiring it into `build_index.R` means
  adding a `group` column to the gene table and keying search on that instead.
- One canonical TSS per gene. Alternative promoters are collapsed. The
  transcript-level TSSs are computed during the build but not yet exposed —
  adding an isoform selector means surfacing `tss_all` from the gene table.
- Brush-zoom sets the window from the dragged span and re-centres pan. It does
  not undo, beyond double-click resetting to ±5 kb.
- No peak track. A narrowPeak BED per species, tabix indexed the same way as
  the models, would slot into the band below the signal with about 20 lines.

## Untested

These scripts have not been run — there is no R in the environment they were
written in. Check first: that your GTFs expose `gene_name` (Ensembl and
GENCODE do; some NCBI GTFs use `gene` instead, which the parser falls back to),
and that `build_index.R` reports sensible sequence-name matches for every
species before you start the app.
