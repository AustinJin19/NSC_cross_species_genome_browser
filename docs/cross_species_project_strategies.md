# Cross-species regulatory genomics: project strategies

Draft date: 2026-10-09

Status: Planning document covering both proposed strategies; no model training or browser changes are performed by this document.

The two strategies expand different dimensions of the project: **strategy 1 adds predicted genomic modalities to species with measured ATAC; strategy 2 adds predicted accessibility for mammals without measured ATAC.**

## Project foundation

The core resource is ATAC-seq across 35 species, paired with their genome sequences. The proposed analysis should begin with a comparable tissue, such as liver, and preserve biological replicates. Species-specific gene annotations and gene orthology assignments are also required for gene-level comparisons.

The objective is to expand the functional interpretation of accessible regions: which genes they may regulate, what regulatory activity they may support, and how those relationships differ across species.

### Shared comparative framework — HALPER + IPP

**HALPER + IPP remains part of the project for finding corresponding regulatory regions across species and assessing conservation.** It supports both strategies below. The detailed implementation remains in [Combining HALPER and IPP for regulatory-element discovery in 35 species](HALPER_IPP_regulatory_workflow.md).

Use the complementary mappings to construct candidate orthologous regulatory-region groups from the multispecies ATAC catalog. Retain agreement, tool-specific candidates, ambiguous mappings, and mapping failures separately. A union of candidates can improve recovery, but each candidate still requires mapping-quality assessment; tool agreement alone does not establish biological conservation.

Evaluate two separate properties of each group:

- **Genomic correspondence:** Is there credible evidence that the regions are orthologous? Record alignment, positional, and mapping-quality evidence.
- **Accessibility conservation:** Are the corresponding regions accessible in the same tissue across species? Use measured ATAC in the 35 profiled species and explicitly labeled predictions in additional mammals.

This produces a shared catalog distinguishing orthologous regions with conserved accessibility, orthologous regions with divergent accessibility, candidate lineage-restricted elements, and unresolved cases. Failure to map is not sufficient evidence for a lineage-specific gain or loss. Accessibility conservation also does not by itself prove conservation of the target gene or enhancer function.

For **strategy 1**, this catalog enables comparison of candidate enhancer–gene links at corresponding elements. For **strategy 2**, it defines comparable loci for sequence-based prediction and helps group homologous examples during training/test splitting. The common deliverable is a region-by-species table containing coordinates, mapping provenance/confidence, measured or predicted accessibility, and missing-data status.

## Strategy 1 — Regulatory networks of genes across species

### Rationale and proposed aim

Reconstruct putative enhancer–gene regulatory networks across species by combining measured ATAC-seq accessibility with DNA sequence and computational predictions of complementary genomic modalities. This approach can expand the breadth of epigenetic information available for species with limited experimental resources, enabling prioritization of candidate regulatory changes for follow-up.

EPCOTv2—the software referred to as EPVCOTV2 in the initial proposal—uses DNA sequence and ATAC-seq to predict transcriptional, epigenomic, transcription-factor-binding, and three-dimensional chromatin signals. These outputs provide several ways to investigate accessible elements beyond their accessibility alone. [EPCOTv2 repository](https://github.com/liu-bioinfo-lab/general_AI_model)

This increases **computational modality coverage**, rather than physical sequencing depth. Predicted RNA signals are not measured RNA-seq or automatically equivalent to TPM. Predicted contacts do not by themselves establish functional enhancer–gene regulation, and the model does not directly provide protein abundance. The reported human/mouse work does not establish accuracy in all 35 species; transfer must be evaluated. [EPCOTv2 paper](https://academic.oup.com/nar/article/53/21/gkaf1269/8340987)

### Questions to address

- Which accessible distal regions are plausible regulators of each gene?
- Do orthologous genes retain similar regulatory inputs across species, or show changes in their candidate enhancer repertoire?
- Which regulatory differences are reproducible across biological replicates and robust to alternative prediction methods?
- If relevant to the project, are gene-level regulatory features associated with maximum lifespan after accounting for phylogeny and technical covariates?

### Required inputs

| Input | Purpose |
|---|---|
| Genome FASTA and chromosome sizes for each species | Sequence input and coordinate validation |
| ATAC alignments or suitably processed quantitative signal, plus peaks | Model input and candidate regulatory regions; peak BED files alone do not contain the full quantitative input |
| Gene/transcript annotation on the same assembly | Promoters, alternative transcription start sites, and target gene IDs |
| Gene orthology table | Comparable gene identities across species |
| Tissue, replicate, sex, age, protocol, and assembly metadata | Quality control and interpretation of biological versus technical differences |
| Compatible pretrained model and recorded preprocessing configuration | Reproducible inference |
| Independent measured RNA, histone, contact, or perturbation data where available | Evaluation; not assumed to exist for every species |

Pretrained inference and model training are different activities. Running a compatible checkpoint does not require generating all target assays first. Training or fine-tuning against those assays requires suitable measured labels.

### Proposed workflow

1. **Build a consistent input manifest.** Match all files to the exact assembly and tissue. Assess ATAC quality, depth, mappability, and replicate agreement. Use the checkpoint's expected signal transformation, with appropriate species-specific genome parameters, rather than copying human constants into every species.

2. **Define candidate regulatory elements and promoters.** Construct reproducible accessible regions within each species. Keep promoter-proximal and distal candidates distinguishable. Retain alternative TSS annotations where supported; assigning every peak to the nearest gene is a baseline, not the final network.

3. **Pilot multimodal predictions.** Begin with a manageable set of loci in human, mouse, and NMR. Record checkpoint, inference window, output resolution, and normalization. Store predictions separately from measured tracks. Evaluate edge effects and consistency between overlapping inference windows before genome-wide processing.

4. **Generate candidate enhancer–gene edges.** Pair accessible distal elements with annotated promoters within the method's supported range. Summarize predicted contact at their anchors and retain distance, accessibility, and promoter identity. Label these as contact-supported candidate links, rather than confirmed regulatory interactions.

5. **Prioritize functional effects.** Where supported, perturb candidate sequence and/or accessibility in silico and evaluate the change in the target transcriptional prediction. Specify exactly which input was perturbed. Use matched control regions and test sensitivity to perturbation size; large artificial perturbations can produce unreliable out-of-distribution inputs.

6. **Compare with alternative link predictors.** Use an accessibility-and-distance baseline, an appropriate ABC configuration, and a direct link predictor where its required reference resources are available. Keep method-specific scores instead of treating their numerical scales as interchangeable probabilities. If predicted contacts or histone signals are used inside an ABC-like calculation, label it as an adapted model and calibrate its transformations and thresholds.

7. **Build species-specific networks.** Represent candidate regulatory regions and genes as nodes, with scored region-to-gene edges. Optionally add TF-to-region motif evidence, but distinguish a motif match from measured TF occupancy. Preserve replicate evidence and uncertain or missing assignments.

8. **Compare networks through gene and regulatory-region orthology.** Compare regulatory inputs to one-to-one orthologous genes first, and use the shared HALPER + IPP catalog to compare links at corresponding regulatory elements. Multiple-copy genes and ambiguous region mappings require explicit handling. Similar target genes do not establish enhancer homology. A gene-centered comparison can still assess regulatory architecture when the contributing enhancers cannot be mapped confidently.

### Other software and suitability for the available data

The methods below have different input requirements. There is no assumption that every listed tool directly accepts bulk ATAC plus FASTA and produces validated links in any species.

| Tool | Inputs and additional requirements | Output and role in this project |
|---|---|---|
| **EPCOTv2** | DNA sequence plus quantitative ATAC input and a compatible checkpoint | Multimodal signal predictions; derive candidate links using promoter contacts and perturbation analysis. Primary proposed prediction engine. [Source](https://github.com/liu-bioinfo-lab/general_AI_model) |
| **ENCODE-rE2G** | Accessibility, gene/candidate annotations, and contact/reference features; model options depend on ATAC versus DNase and optional H3K27ac | Direct enhancer–gene scores. Particularly relevant to an ATAC-led analysis, but the released workflow recommends pooled megamap Hi-C when cell-specific Hi-C is unavailable and does not provide trained power-law models. Human reference resources cannot simply be reused in NMR coordinates. [Source](https://github.com/EngreitzLab/ENCODE_rE2G) |
| **ABC model** | Accessible candidates, activity estimates, gene annotations, and contact estimates; exact assay requirements depend on configuration | Enhancer–gene scores based on activity × contact. An ATAC-led configuration with a distance-based contact approximation provides an exploratory baseline; the approximation is not measured 3D structure. [Methods](https://abc-enhancer-gene-prediction.readthedocs.io/en/stable/usage/methods.html) |
| **Enformer** | DNA sequence and pretrained model; ATAC can select candidates but is not the model's input | Expression/epigenomic predictions and sequence perturbation or attribution scores. Useful as a sequence-only comparator, not a direct readout of the supplied sample's ATAC state. [Paper](https://www.nature.com/articles/s41592-021-01252-x) |
| **SPEID** | Candidate enhancer and promoter sequences, with a pretrained or appropriately trained interaction classifier | Sequence-based enhancer–promoter interaction predictions. ATAC can define candidates upstream. A legacy comparator whose software compatibility and tissue/species transfer require checking. [Repository](https://github.com/ma-compbio/SPEID) |
| **Cicero** | Single-cell ATAC across many cells | Coaccessibility links that can connect distal sites with promoters. Not appropriate for a few bulk ATAC replicates; included for comparison, not as a proposed R-based implementation. [Project](https://cole-trapnell-lab.github.io/projects/cicero/) |
| **SCARlink** | Matched single-cell ATAC and RNA measurements | Gene-specific regulatory-region importance and expression predictions. Requires additional multiome data and is not an ATAC-only option for the current bulk dataset. [Repository](https://github.com/snehamitra/SCARlink) |

**Recommended starting combination:** EPCOTv2 for multimodal hypotheses, a distance/accessibility baseline for context, and ABC or ENCODE-rE2G where their inputs and calibration are defensible. Keep Enformer as an optional sequence-only comparison. Agreement among methods is supporting evidence, but shared inputs and training data mean it is not independent experimental validation.

### Evaluation and safeguards for comparative inference

- **Validate against independent measurements.** Compare predicted RNA with measured RNA and predicted contacts with measured contacts where tissue and assembly match. Comparing two outputs of the same model is a consistency analysis, not independent validation.
- **Prevent leakage.** Hold out chromosomes or genomic regions during fitting and calibration. For cross-species evaluation, also consider held-out species/clades and homologous regions shared with training data.
- **Measure useful performance.** For links with experimental labels, report precision–recall and performance stratified by enhancer–promoter distance. For quantitative assays, assess rank agreement and calibration. Compare with simple baselines.
- **Separate missing evidence from biological absence.** Assembly gaps, annotation gaps, low coverage, and failed model transfer must not be called lineage-specific regulatory losses.
- **Control comparability.** Do not compare raw predicted scores across species without checking normalization and calibration. Tissue cell composition and sequencing protocol can imitate evolutionary changes.
- **Keep tissue provenance explicit.** The existing NMR hypothalamus ABC test data cannot validate liver enhancer–gene predictions. Existing model-predicted Micro-C loops are also predictions, unless separately established as experimental measurements.
- **Use phylogenetic analysis for trait associations.** If testing lifespan relationships, account for phylogeny, multiple testing, technical covariates, and uncertainty in the predicted regulatory features.

### Deliverables and staged execution

1. **Pilot:** input manifest, quality-control report, representative locus predictions, and baseline comparisons. Identify where species transfer is unsupported before scaling.
2. **Network release:** a per-species edge table with species, assembly, tissue, sample, enhancer coordinates, target gene/TSS, distance, method, score, evidence type, and model version. Store measured, predicted, and independently validated evidence separately.
3. **Comparative analysis:** ortholog-centered summaries of candidate enhancer number, accessibility, and link support, together with missingness and uncertainty. Add element-level conservation only where orthology is established.
4. **Browser presentation:** optional predicted expression tracks and candidate regulatory arcs with visible tissue, model, and evidence labels. Preserve the underlying ATAC measurements.
5. **Follow-up shortlist:** robust candidate regulatory changes for independent assays or perturbation experiments.

## Strategy 2 — Expand mammalian coverage with sequence-based accessibility prediction

### Rationale and proposed aim

Use the ATAC-seq data from 35 species as training labels for a machine-learning model that predicts tissue-specific chromatin accessibility from DNA sequence. Apply the trained model to additional mammalian genomes to expand the phylogenetic breadth of the study beyond species with experimental ATAC-seq.

The closest precedent to the proposed Zoonomia comparison is **TACIT: Tissue-Aware Conservation Inference Toolkit**. TACIT learned tissue-specific open-chromatin predictions from sequence and associated predictions across mammals with phenotypes, including brain-size residuals. Its principal result was enhancer–phenotype association, rather than simply predicting brain size directly from a genome. Our proposed adaptation is to learn the liver accessibility code from the 35-species resource and investigate relevant mammalian traits. [TACIT paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10322212/)

TACIT's repository documents sequence-model training and subsequent phylogenetic association analyses, including links to brain and liver models. It is a methodological starting point; the proposed 35-species model still needs its own training, validation, and tissue-specific interpretation. [TACIT software](https://github.com/pfenninglab/TACIT)

### What the model would learn

During training, the input is a genomic DNA sequence and the target is accessibility measured by ATAC in the corresponding species and tissue. During prediction in an additional species, the input is its DNA sequence; experimental ATAC from that target species is not required.

Start with a clear, achievable target: **the probability that a candidate region is accessible in liver**. Quantitative normalized accessibility or base-resolution signal can be a later objective if assay comparability supports it. Neither an accessibility probability nor a sequence-based signal prediction should be presented as newly measured ATAC read counts.

With bulk liver training data, the output reflects the sampled bulk tissue. It cannot resolve unmeasured cell types or reliably infer effects of age, environment, or treatment that are not encoded in sequence and represented in the training design.

### Proposed workflow

1. **Harmonize the measured training resource.** Process the 35 species consistently, documenting tissue, replicate quality, coverage, assembly quality, and mappability. Define reproducible positives and distinguish reliably inaccessible regions from regions with insufficient evidence.

2. **Construct sequence examples and informative negatives.** Extract fixed-length windows around candidate regulatory regions. Include accessible elements and suitable confidently inactive examples, balancing GC content, repeat content, genomic context, and species representation. Avoid teaching the model only promoter proximity or peak-calling artifacts. Treat ambiguous labels as uncertain rather than forcing them into the negative class.

3. **Establish orthology using HALPER + IPP.** Apply the shared comparative workflow to define candidate corresponding regulatory regions and track mapping quality. Use these groups to assess conserved versus divergent accessibility and to compare sequence-based predictions across species. Preserve unalignable, ambiguous, or duplicated regions explicitly; do not treat them as inaccessible by default.

4. **Train a sequence-only model.** Begin with a pooled multispecies convolutional neural network and a simple sequence-feature baseline. Evaluate species-balanced sampling so well-sampled species do not dominate. More complex architectures or quantitative prediction should be justified by improved held-out performance, not chosen solely for model size.

5. **Test transfer before extrapolating.** Evaluate held-out species and entire clades. Also group homologous loci when splitting genomic training and test examples. Report two distinct settings: prediction at orthologs of loci represented during training, and prediction at previously unseen locus families. A strict joint species-and-locus holdout tests the harder form of generalization. No test-species ATAC should enter fitting, threshold selection, or calibration.

6. **Predict accessibility in additional genomes.** First score orthologs of a multispecies union of candidate elements. This yields an interpretable region-by-species matrix. A subsequent genome-wide scan can seek new candidates absent from that union, but needs separate false-positive assessment. A reference-only candidate list cannot discover every lineage-specific element.

7. **Associate regulatory predictions with traits.** For traits such as maximum lifespan, test whether variation in predicted accessibility is associated with the trait using phylogeny-aware models. Consider body mass and other justified covariates, correct for multiple testing, and assess robustness across clades. Use measured accessibility in the original species as an empirical check where possible.

8. **Prioritize and validate candidates.** Identify robust candidate elements, investigate nearby or independently linked genes, and select representative species or loci for experimental ATAC, reporter assays, or other appropriate validation. An association with a trait is not proof that the element causes the trait difference.

### Evaluation and interpretation

| Question | Suggested assessment |
|---|---|
| Does the model recognize accessible regions? | Precision–recall, calibration, and recall at a controlled false-positive rate; report class prevalence and negative sampling |
| Does it transfer across evolution? | Leave-species-out and leave-clade-out performance, stratified by distance to training species |
| Does it predict regulatory change rather than only sequence conservation? | Accuracy on orthologous loci with measured accessibility differences; comparison with conservation and sequence-similarity baselines |
| Are predictions stable? | Replicate-label sensitivity, model ensembles, and performance across GC, repeats, promoter/distal status, and assembly quality |
| Are trait associations robust? | Phylogenetic null analyses, multiple-testing correction, clade sensitivity, and checks for assembly/batch confounding |

Thirty-five species provide many genomic training examples, but those examples are correlated through shared ancestry and homologous loci. Likewise, hundreds of predicted species are not hundreds of independent experimental ATAC datasets. Prediction uncertainty and phylogenetic dependence must remain part of the interpretation.

An absent prediction is not automatically an inaccessible region: it may reflect missing sequence, failed orthology mapping, poor model transfer, or a low-confidence score. Keep these categories distinct from predicted biological inactivity.

### Deliverables

- A reproducible training dataset with sequence, accessibility labels, species/tissue metadata, and split assignments.
- A validated sequence-to-accessibility model with documented transfer limits and uncertainty estimates.
- A matrix of candidate regulatory regions × mammalian species containing predicted accessibility, mapping status, and confidence.
- Phylogeny-aware regulatory-element–trait association results, if trait analysis is pursued.
- Browser tracks clearly labeled as predicted accessibility, with measured tracks retained separately.

## How the two strategies connect

HALPER + IPP supplies the shared regulatory-region correspondence and conservation framework for both strategies. Multimodal prediction adds functional hypotheses to this catalog, while sequence-based accessibility prediction extends its coverage to additional mammals.

| Dimension | Strategy 1: expand modalities | Strategy 2: expand species |
|---|---|---|
| Main input | Measured ATAC plus DNA sequence | DNA sequence at inference; 35-species ATAC supplies training labels |
| Main output | Predicted complementary signals and candidate enhancer–gene networks | Predicted tissue-specific accessibility in additional mammals |
| Immediate scope | Species with suitable measured ATAC | Additional species with suitable genome sequences |
| Main uncertainty | Multimodal model transfer and functional interpretation of links | Sequence-model transfer and unobserved tissue/state effects |

The strongest initial integration is to connect trait-associated elements from strategy 2 to gene hypotheses supported by strategy 1 in experimentally profiled species. Transfer gene assignments to other species only with appropriate orthology and regulatory evidence.

Feeding predicted ATAC from strategy 2 into EPCOTv2 is a possible later experiment, not an automatically valid extension. It compounds prediction error and may create an input distribution unlike measured ATAC. It would require a held-out benchmark comparing this cascade against the measured-ATAC pipeline before biological conclusions are drawn.

Together, the strategies aim to increase the **functional depth and phylogenetic breadth** of the project while retaining a clear distinction between measurements, predictions, and validated regulatory mechanisms.
