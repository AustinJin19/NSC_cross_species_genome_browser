# Project folder guide

The local dataset was reorganized by data type on 2026-10-10. Files were moved, not copied or transformed. The browser and Python preprocessing/analysis scripts use the new locations.

| Folder | Contents |
|---|---|
| `reference_genomes/` | FASTA genomes, matching `.fai` indexes, GTF annotations and `.tbi` indexes, grouped by species; hg19 chromosome sizes |
| `data/atac/` | ATAC bigWigs and matching peak calls, grouped by species |
| `data/atac/mouse_legacy/` | Preserved older mouse assembly data; the browser continues to use `mouse_mm10/` |
| `data/atac/mole_rat_peaks/` | NMR, BMR, and DMR samples; NMR1 is preserved on disk but excluded from the browser |
| `data/chipseq/mhetglav3/` | NMR H3K4me3 and H3K27ac regions |
| `data/chipseq/ChIP_seq/` | Public histone peaks and their self-contained liftover pipeline, chains, results, and provenance |
| `data/interactions/human_hic/` | Human hg38 `.mcool` matrix (active); original `.hic` retained as an archive |
| `data/interactions/hicorr/` | Original HiCorr BEDPE interactions |
| `data/interactions/nmr_microc/` | NMR predicted Micro-C BEDPE files and metadata, including enrich1.0 |
| `data/interactions/nmr_abc/` | Hypothalamus ABC test BEDPE input |
| `data/regulatory_annotations/ENCODE4/` | Human cCRE registry and metadata |
| `data/regulatory_annotations/repeats/` | Original NMR TE repeat BED |
| `data/regulatory_annotations/pgls/` | Original hg19 PGLS LAR BED |
| `data/omics/rna/` | CrossSpecies_Count.xlsx and rnaLongTPM.xlsx |
| `data/omics/protein/` | protLongDF.xlsx |
| `data/omics/lifespan/` | AnAge archive |
| `data/metadata/mole_rats/` | Original species XML metadata |
| `archives/` | Original ChIP download archive and MLS pipeline archive |
| `cross_species_visual_tool/web_app/` | Active Python browser source |
| `cross_species_visual_tool/index/` | Generated transcript indexes and human anchors |
| `cross_species_visual_tool/web_app/data/` | Generated browser SQLite, lifespan JSON, and indexed TE data |
| `HiCorr_pgls_LAR_intersection/` | HiCorr/LAR analysis scripts and results |
| `MSA_transformer/` | Existing analysis scripts, input table, and outputs, kept together |
| `output/` | Exported PDFs, source archive snapshot, migration manifest, and legacy figures |
| `docs/` | Project documentation and research plans |
| `tmp/` | Existing temporary working files |

Reference subfolders retain familiar species-folder names: `selected3` is human hg38, `mouse_mm10` is the current mouse genome, and `mouse_legacy` is the old mouse reference. The three mole-rat references use their scientific names. All original filenames are preserved.

## Verification and provenance

`output/file_organization_manifest.json` records every moved file's old path, new path, byte size, and modification time. All 233 files retained their byte sizes and modification times during the move. Afterward, `data/chipseq/ChIP_seq/liftover/run_liftover.py` was edited to use the new reference paths; the other 232 moved files remain unchanged. Moves used filesystem renames rather than copying large genomes or the Hi-C matrix.

Historical result metadata and archived source ZIPs are preserved as snapshots and may mention old paths. The legacy R `config.tsv` now points to the reorganized inputs; existing RDS caches may still contain historical absolute paths and require rebuilding if the old R app is used. The active Python browser does not use those RDS caches.

## Start the browser

From the project root:

```bash
python3 cross_species_visual_tool/web_app/server.py
```

Then open http://127.0.0.1:8765/. No data reprocessing is required for this reorganization. Generated indexes stay in their existing application folders.
