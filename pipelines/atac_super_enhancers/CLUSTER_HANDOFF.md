# Continue this analysis with Codex on the cluster

Suggested task message:

> Run the pipeline in `pipelines/atac_super_enhancers` to identify ATAC-defined
> candidate super-enhancers in my species. Read its README first. Locate the
> project root and the filtered BAM directory I provide. Complete and verify the
> BAM mapping in `samples.cluster.tsv`; never choose an ambiguous match or change
> assemblies. Keep sample-level calls separate. Start with the pinned ROSE
> defaults (12.5 kb stitching, 2.5 kb TSS exclusion), test one real sample in a
> compute allocation, inspect its rank curve and validation report, then submit
> the remaining samples using the available scheduler. Use matching-assembly
> exclusion BEDs only when available; document missing masks and QC. Summarize
> per-sample counts and regions supported in at least two samples per species.
> Preserve source data and the current website. Report failures and ambiguous
> inputs instead of silently dropping samples. The local work prepared and tested
> the pipeline; no real BAM scoring has been performed yet.

Required information to provide there:

- Cluster project directory (containing `data/`, `cross_species_visual_tool/`, `pipelines/`).
- Filtered BAM directory and any sample-to-filename naming differences.
- Scheduler partition/account or an existing compute allocation.
- Optional species-specific blacklist/low-mappability BED paths.

The current manifest includes 89 samples / 34 species. Human, macaque and bowhead
whale are included for SE identification; previous MLS plot exclusions do not
apply automatically. Low-TSS samples are flagged, not automatically dropped.
