# Liver ATAC bundle, imported 2026-10-15

Data are grouped by species within the existing project categories:

- `reference_genomes/liver_atac_bundle/<Scientific_name>/`: compressed FASTA,
  FASTA `.fai`/`.gzi`, GTF, and GTF tabix indexes.
- `data/atac/liver_atac_bundle/<Scientific_name>/`: sample CPM bigWigs and narrowPeak files.
- `data/liver_atac_bundle/`: original source manifest and validation report,
  portable `manifest.local.tsv`, import report and browser validation report.

The source manifest records upstream cluster paths and is retained unchanged.
Use `manifest.local.tsv` for paths relative to the project root. Files were moved,
not expanded or duplicated. The `.fa.gz` references are indexed compressed FASTA
and can be read with pysam directly.

The importer adds species with signal and annotation files. Existing browser
species retain their established tracks, models, auxiliary annotations and
activity scores. After the 2026-10-16 cleanup, identical bundle copies resolve to
the retained original tracks through `data/liver_atac_bundle/path_overrides.json`.
The two corrected ASM peak files remain in the bundle directory; the incompatible
old ASM peaks were removed. NMR1 signal and peak files were removed because that
sample is excluded; its local manifest entries are `-`. The source manifest and
original validation reports remain historical records.
Four species have no supplied bigWig and are not added:
Cricetomys gambianus, Scalopus aquaticus, Ondatra zibethicus, Tamias striatus.
Their reference and peak files are retained.

Copy `path_overrides.json` with the other bundle metadata when deploying. The
importer honors these mappings so re-importing does not replace retained paths
with missing duplicate paths or reintroduce the excluded NMR1 sample. The cleanup
record is `output/data_cleanup_2026-10-16.json`. Older mouse ATAC tracks are retained
because the legacy R browser configuration still references that assembly.

New species are configured in `web_app/data/additional_species.json` and use
portable `index/<Scientific_name>_models.tsv.gz` plus tabix indexes. The server
loads these additions when present; installations without the bundle retain their
original catalog. Restart the server and refresh the page after import.

Rebuild with the Python environment containing pyBigWig and pysam:

```bash
python cross_species_visual_tool/web_app/import_liver_bundle.py
```

Local bundled dependencies, if used, require
`PYTHONPATH=cross_species_visual_tool/web_app/.deps`.

Validation checks every new bigWig chromosome name and length against the supplied
FASTA index, peak bounds against that index, and annotation bounds against signal
chromosomes. A FASTA fetch validates compressed random access. This checks internal
coordinate compatibility; unspecified assembly accessions remain unspecified.

Symbols come from explicit gene_name/gene attributes, or the source bed2gtf
transcript-ID symbol field. Bowhead ORF annotations use annotation_label with the
trailing source-species suffix removed (e.g. HAS2_CISTUGOS becomes HAS2); full
source labels remain in the transcript ID and original GTF. These homology-derived
labels and protein mnemonics are not validated one-to-one orthologs.

Gene activity now includes new species via an indexed, on-demand calculation.
Run `python cross_species_visual_tool/web_app/activity.py --register-additions`
after importing, followed by `python cross_species_visual_tool/web_app/build_mls.py`.
All six activity strategies are supported. Lifespan metadata use the existing local
references; unmatched species are excluded from correlations. The liver RNA/protein
section and earlier offline correlation/enrichment reports remain unchanged.
