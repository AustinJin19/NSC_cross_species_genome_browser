# ATAC-defined candidate super-enhancers across species

This pipeline runs on the VM/HPC using **filtered BAMs**, the existing per-sample
ATAC peaks, and the exact species annotation indexes used by the genome browser.
It does not use CPM bigWigs for scoring. The shipped template covers **89 samples
from 34 species**, including human, macaque and bowhead whale. These are SE calls,
not MLS correlations; correlation-specific species exclusions are not inherited.
Subset the manifest if a different species or sample set is wanted.

**No real BAM analysis has been run locally.** BAM paths in `samples.template.tsv`
are deliberately blank. Local peak/annotation checks are recorded in
`local_input_validation.json`. The workflow does not change the website or source data.

## Cluster quick start

Copy this entire folder to `pipelines/atac_super_enhancers/` under the cluster's
project root. The cluster also needs the referenced `data/atac/` peak files and
`cross_species_visual_tool/index/*_models.tsv.gz` files in the same relative layout.
`reference_sizes/` is included; the bigWigs and FASTAs are not needed to run the
shipped manifest. All non-absolute input paths are relative to `--project-root`.

```bash
cd /path/to/CS_Genome_Browser
conda env create -f pipelines/atac_super_enhancers/environment.yml
conda activate atac-rose

python pipelines/atac_super_enhancers/workflow.py discover-bams \
  --manifest pipelines/atac_super_enhancers/samples.template.tsv \
  --bam-root /path/to/filtered_bams \
  --output pipelines/atac_super_enhancers/samples.cluster.tsv
```

**Review `samples.cluster.tsv` before submitting.** Discovery only matches a
sample-name prefix followed by a delimiter. It cannot determine whether a file
is the correct filtered BAM or assembly. Zero/multiple matches stay blank and
are listed in `.unresolved.json`. Edit those paths explicitly. Different sample
naming conventions require manual mapping. Never select one ambiguous BAM arbitrarily.

First test the software and run one real sample on a compute node:

```bash
python -m unittest discover -s pipelines/atac_super_enhancers -p test_workflow.py
python pipelines/atac_super_enhancers/smoke_test.py

python pipelines/atac_super_enhancers/workflow.py run-sample \
  --manifest pipelines/atac_super_enhancers/samples.cluster.tsv \
  --project-root "$PWD" \
  --output outputs/atac_rose_12p5kb \
  --sample NMR3Liver
```

Inspect that sample's `run.log`, `bam_validation.json`, `summary.json` and
`ATAC_Plot_points.png` before scaling up. Do not run all samples on a login node.

For SLURM, activate the environment before submission:

```bash
bash pipelines/atac_super_enhancers/submit_slurm.sh \
  pipelines/atac_super_enhancers/samples.cluster.tsv \
  "$PWD" outputs/atac_rose_12p5kb
```

The array defaults to four concurrent jobs, one CPU, 24 GB RAM and 24 hours per
sample. Edit `run_array.sbatch` for the site's partition/account/time/memory rules;
adjust `MAX_PARALLEL_JOBS` for concurrency. The scripts assume the activated conda
environment is inherited by the compute nodes. Use an interactive allocation to
check this on your cluster. For a non-SLURM allocation, `run_local.sh` accepts the
same three arguments and runs samples sequentially.

After all jobs finish:

```bash
python pipelines/atac_super_enhancers/workflow.py summarize \
  --manifest pipelines/atac_super_enhancers/samples.cluster.tsv \
  --output outputs/atac_rose_12p5kb
```

A completed sample is skipped only when the inputs, settings, pipeline version,
and output checksums match. Incomplete/changed output directories are not
overwritten; inspect the failure and use a new output directory for a retry.
Summary refuses incomplete manifests rather than silently dropping failed samples.

## Analysis definition

| Step | Default |
|---|---|
| Unit of analysis | Each sample independently, in its own native assembly |
| Constituent regions | Existing sample narrowPeak/BED calls, exact duplicate intervals removed |
| Stitching distance | 12,500 bp; configurable with `--stitch` |
| TSS exclusion | ROSE `-t 2500` behavior, using all indexed transcript TSSs |
| Ranking signal | Pinned ROSE BAM mapper: reads extended by 200 bp, RPM normalization, floor=1 |
| SE threshold | Unmodified upstream `ROSE_callSuper.R` tangent/rank-curve cutoff |
| Candidate label | Score strictly greater than the ROSE cutoff |
| Per-species support | Exact genomic bases called SE in at least two distinct samples |

The TSS rule removes peaks **contained within** an annotated ±2.5 kb TSS region;
it is not an overlap-based promoter subtraction. Upstream ROSE also breaks
stitched regions overlapping more than two distinct gene TSSs (its ±50 bp test)
back into constituent loci. All transcript isoforms are retained; this does not
use the browser's single chosen gene locus or human MANE anchor.

This uses standard ROSE **read-based scoring**, including both alignment mates
and its legacy 200 bp extension/floor, not a new fragment-count/Tn5-insertion
implementation. CPM bigWigs only supplied the reference sequence dictionaries
when creating the manifest. Document these choices in a methods section.

### Assembly and BAM checks

- Coordinate-sorted, indexed paired-end BAMs are required. `.bai` and `.csi` are accepted.
- Every sequence in the supplied reference dictionary must be present in the BAM
  with the same length. Peak and model coordinates must be within those bounds.
  Names and lengths are necessary checks, not proof that two assemblies are identical.
- The first 100,000 alignments must have MAPQ ≥30, be proper pairs and lack
  unmapped, secondary, supplementary, duplicate and QC-failed flags. This is a
  **preflight spot check**, not whole-BAM QC or a filtering step. Supply BAMs that
  were fully filtered upstream. `--min-mapq` changes the required threshold only.
- No BAM is reheadered, rewritten, re-aligned, pooled or automatically filtered.
- `exclude_bed` optionally specifies a matching-assembly BED of low-mappability,
  blacklist or other excluded regions. Any overlapping candidate peak is removed.
  Blank means **no such mask was applied**, which is recorded in the output.
- `exclude_chroms` is a semicolon-delimited list. The template lists conventional
  mitochondrial names and `chrY`; ROSE itself excludes literal `chrY`. For
  accession-named assemblies, add the **verified** mitochondrial/Y accessions
  before running if those sequences should be excluded. No accession identity is guessed.
- Samples with TSS enrichment below 4 are flagged and retained by default,
  following the current all-sample browser policy. QC is unavailable for some samples.
  For a quality-filtered analysis, explicitly subset the manifest; optionally use
  `--min-tss-qc 4` to enforce the filter (unknown QC then fails).

### Version and compatibility

ROSE is vendored from <https://github.com/stjude/ROSE>, commit
`d845b98e8e8fa657b3d3e9af91653a8eea7449b4`, with its license and SHA256 checksums
in `vendor/source.json`. No download is needed on compute nodes.

The runner calls the pinned ROSE Python stitching/scoring functions and the
unchanged R cutoff script synchronously. This avoids the upstream launcher's
background-process/file-existence race and supports failures as nonzero exits.
Two runtime compatibility adjustments preserve native contig names (disable
automatic `chr` insertion/removal) and clamp BAM query starts to 1 at contig edges.
The vendored source files remain unmodified and are checked before execution.

Custom UCSC-layout annotation tables are generated from browser transcript
indexes. Synthetic unique transcript IDs avoid collisions; original IDs are saved
in `annotation_refseq.ucsc.transcript_ids.tsv`. Coding regions are not inferred.
Our downstream annotation reports the gene with the closest TSS to the region
midpoint and all genes with a TSS within 50 kb of its edges. This replaces the
upstream gene-mapper and does **not** establish enhancer–gene contacts or orthology.

ROSE uses legacy coordinate handling internally. Input and final BED files retain
the original BED 0-based, half-open boundaries; its intermediate GFF is the
upstream BED-to-GFF convention, not a generic 1-based annotation GFF. Do not add
an extra coordinate conversion to these intermediate files.

## Outputs

Each `<output>/<species>/<sample>/` contains:

- `candidate_se.bed`: candidate SE regions, BED4.
- `candidate_se.tsv`: candidate coordinates, RPM-integrated signal, rank,
  constituent-peak count, nearest gene and nearby TSS genes.
- `all_regions.tsv`: all ranked regions including typical enhancers.
- `ATAC_Plot_points.png`: ROSE rank plot and threshold.
- Original ROSE tables, mapped signal files and stitching debug output.
- `inputs.json`, `runtime_versions.json`, `bam_validation.json`, `summary.json`,
  `run.log`, `SUCCESS.json`. Inputs other than the BAM are SHA256-checked; large BAMs
  are tracked by path, size and modification time rather than a full-file checksum.

`species_summary/` contains sample/species count tables and, for each species:

- `*_candidate_se_union.bed`: union of any sample's candidate SEs.
- `*_candidate_se_shared_ge2.bed`: bases supported by at least two samples.
- `*_support_segments.tsv`: exact support count and sample identities per segment.

Shared regions are intersected at base resolution, then adjoining qualifying
segments are merged. Transitive chains of overlapping calls are **not** treated
as uniformly replicated. A one-sample species has no replicated calls; this is
recorded as `single_sample_only`. Replicate support is not a fresh ROSE call on
a pooled species BAM and may cover only part of a sample's SE interval.

## Sensitivity and interpretation

For a separate 5 kb stitching analysis, append `--stitch 5000` to a run command
and use a new output directory, e.g. `outputs/atac_rose_5kb`. Do not combine
different stitching settings in one species summary. Compare replicate
consistency and the rank curves before choosing a final analysis.

Call the results **ATAC-defined candidate super-enhancers**. ATAC accessibility
alone does not establish enhancer function or liver specificity. Use H3K27ac or
functional validation where available; human/mouse dbSUPER liver BEDs provide an
independent comparison. Scores and counts across assemblies are affected by peak
calling, depth, quality, annotation completeness and assembly fragmentation.
RPM is not cross-species batch correction. A raw SE count difference is not
evidence of a biological species difference. Orthology/synteny mapping and matched
QC are separate requirements for cross-species comparisons or MLS association.

## References

- ROSE code: <https://github.com/stjude/ROSE>
- Whyte et al., Cell (2013): <https://doi.org/10.1016/j.cell.2013.03.035>
- Lovén et al., Cell (2013): <https://doi.org/10.1016/j.cell.2013.03.036>
- ATAC precedent, dbscATAC (2025): <https://doi.org/10.1093/bioinformatics/btaf364>
- Alternative 5 kb stitching, Wibisana et al. (2022): <https://doi.org/10.1371/journal.pgen.1010235>
