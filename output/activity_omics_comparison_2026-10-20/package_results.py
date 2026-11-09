#!/usr/bin/env python3
"""Consolidate comparison results, document methods, and bundle deliverables."""
import csv
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'cross_species_visual_tool/web_app').is_dir())
HERE=Path(__file__).resolve().parent
OUT=ROOT/'outputs/activity_omics_comparison_2026-10-20'
LABELS={'tss2kb':'TSS ±2 kb','tss5kb':'TSS ±5 kb','body':'Gene body + 2 kb',
        'tss200kb':'TSS ±200 kb','tss500kb':'TSS ±500 kb','rna':'Transcriptome','protein':'Proteome'}


def read(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))


def write(path,rows,fields):
    with path.open('w') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)


def main():
    summary=json.loads((OUT/'summary.json').read_text())
    primary=[r for r in summary['comparisons'] if r['variant']=='mass_adjusted']
    positive=[];negative=[];opposite=[];fdr=[]
    for s in primary:
        m,assay=s['method_a'],s['cohort']
        rows=read(OUT/'mass_adjusted'/assay/f'{m}_vs_{assay}_all.tsv')
        for r in rows:
            r=dict(method=LABELS[m],assay=LABELS[assay],**r)
            if float(r['activity_p'])<.05 and float(r['expression_p'])<.05:
                (positive if r['direction']=='positive_both' else negative if r['direction']=='negative_both' else opposite).append(r)
            if float(r['activity_q'])<.05 and float(r['expression_q'])<.05:
                fdr.append(r)
    fields=['method','assay','gene','n_species','direction','activity_r','activity_p','activity_q','expression_r','expression_p','expression_q','species']
    for name,rows in [('shared_positive_nominal',positive),('shared_negative_nominal',negative),('opposite_direction_nominal',opposite),('shared_FDR',fdr)]:
        write(OUT/(name+'.tsv'),rows,fields)
    comparison=[]
    for r in primary:
        comparison.append(dict(activity_method=LABELS[r['method_a']],expression_assay=LABELS[r['cohort']],
          genes_tested=r['tested_genes'],Spearman_rho=r['spearman_rho'],same_direction_percent=r['same_direction_percent'],
          positive_both_p05=r['nominal_positive_both'],negative_both_p05=r['nominal_negative_both'],
          opposite_both_p05=r['nominal_opposite_direction'],same_direction_both_FDR05=r['FDR_same_direction'],
          minimum_species=r['minimum_species'],median_species=r['median_species'],maximum_species=r['maximum_species']))
    write(OUT/'gene_overlap_counts.tsv',comparison,list(comparison[0]))
    species=read(OUT/'species_inclusion.tsv')
    n_rna=sum(int(r['RNA_species_present']) for r in species)
    n_protein=sum(int(r['protein_species_present']) for r in species)
    report=['MLS association comparison — 20 October 2026','',
      'Primary analysis: log10 MLS; log10 median ATAC activity; body-mass correction ON.',
      'Gene-body ATAC additionally adjusts for log10 annotated gene length.',
      f'{len(species)} eligible ATAC species; {n_rna} overlap RNA and {n_protein} overlap protein at the dataset level.',
      'Per-gene species counts are smaller when a gene or sufficient genomic region is unavailable.','',
      'Method\tExpression\tGenes tested\tSpearman rho\tShared positive p<.05\tShared negative p<.05\tOpposite p<.05\tShared same-direction FDR<.05']
    for r in comparison:
        report.append(f"{r['activity_method']}\t{r['expression_assay']}\t{r['genes_tested']}\t{r['Spearman_rho']:.3f}\t{r['positive_both_p05']}\t{r['negative_both_p05']}\t{r['opposite_both_p05']}\t{r['same_direction_both_FDR05']}")
    report += ['', 'How to read these results',
      'Spearman rho compares the gene-by-gene MLS correlation coefficients between two assays. It is not a direct correlation of ATAC activity with RNA/protein abundance.',
      'A shared positive or negative gene has p<0.05 in both assays and the same association sign. Nominal matches are exploratory; use the separate FDR columns for multiplicity-adjusted evidence.',
      'Each RNA comparison uses exactly the same genes and per-gene species across RNA and all five ATAC methods. Protein comparisons are matched separately. RNA and protein cohort sizes therefore differ.',
      'The joint correlation map uses one common cohort across all seven methods; its coefficients differ from the two expression-specific comparison columns.',
      'The gene heatmap displays up to 50 genes selected by the smallest worst p value in an ATAC/expression same-direction pair in the joint cohort. Asterisks indicate nominal p<0.05. It is a selected display, not independent validation.','',
      'Sensitivity analyses',
      'without_mass/: body mass omitted from all fits; gene-body length correction retained.',
      'without_ASM/: African spiny mouse omitted because RNA source common/species labels disagree.',
      'at_least_10_species/: require at least 10 matched species per gene.',
      'mass_adjusted/all_available/: per-method reference using its own available species; do not confuse with matched heatmap counts.']
    (OUT/'RESULTS.txt').write_text('\n'.join(report)+'\n')
    methods='''DATA AND METHODS

ATAC: current browser gene_activity.sqlite loci and current species bigWigs. Exact piecewise-constant bigWig integrals; signal outside represented intervals contributes zero, absent chromosomes are missing. TSS ±2/±5/±200/±500 kb means total spans 4/10/400/1000 kb. Gene-body method includes the complete annotated body plus 2 kb upstream, respecting strand. Outer locus TSS is the left body edge on + and last body base on -. Coordinates and clipping match the browser. All annotation biotypes and gene widths are retained.

For each locus, take the median across eligible sample tracks; select the highest-median locus per symbol and species separately for each method. Then exclude the selected locus if fewer than 50% of the requested bases fall inside the assembly chromosome/scaffold bounds in any measured sample. There is no fallback to a lower-activity locus. This availability rule measures assembly span, not read coverage, N content or peak occupancy. Apply log10 to the species median; exclude nonpositive values, with no pseudocount. This analysis uses all signal in the requested regions, not only signal within peaks and not peak counts.

RNA: liver_omics.sqlite, CrossSpecies_Count.xlsx plus the existing rnaLongTPM.xlsx mouse supplement; species median of supplied log.tpm. Protein: protLongDF.xlsx; species median of supplied quant.norm, whose normalization basis is log.quant. These already-transformed source values are not logged again. Samples are aggregated to species; biological samples are not independent observations in the MLS fit. Duplicate assay/species/sample/gene rows were checked and absent. BMR_Turk_Lu is excluded. RNA/protein source workbooks were checked against index provenance and were not stale.

Species: project MLS and mass metadata. Human, macaque and bowhead are excluded throughout, as are missing, ambiguous or lower-bound MLS values. Only species in the current ATAC collection enter expression comparisons. Ellobius lutescens RNA is mapped to Bramus lutescens using the existing project taxonomic mapping. Acomys RNA has a conflicting common-name field (A. russatus) versus species field (A. cahirinus); primary analysis follows the source species field, with a separate ASM-excluded sensitivity. BMR lifespan is the existing Spalax-complex proxy, documented in species_inclusion.tsv.

Primary model: response ~ log10(MLS) + log10(body mass); gene-body ATAC adds log10(gene length in kb). RNA/protein use their supplied transformed abundance as response; ATAC uses log10 median integrated activity. Residualize response and log10 MLS against the covariates, then calculate Pearson partial r and two-sided t p with df = n - number_of_covariates - 2. At least five species are required. Rank-deficient or constant fits are excluded. The sensitivity without mass retains only the gene-length covariate for gene-body ATAC; other methods then use ordinary Pearson correlation. Independent full-model OLS t tests validate the partial-r calculations.

Matching: for each gene, intersect eligible species across the five ATAC methods and RNA (or protein). All methods in that comparison use the exact same species. Retain genes with valid fits in every method. Therefore even a narrow-window comparison is restricted to species with adequate wider-window data; this is intentional to compare window methods fairly. The separate joint analysis intersects all seven datasets. Genes are matched by uppercase symbol only; the current index does not establish one-to-one orthology.

Multiple testing: Benjamini-Hochberg correction across all tested genes separately for each method, matched cohort and sensitivity analysis. A nominal overlap requires both p<0.05; an FDR overlap requires both q<0.05. Positive and negative overlap require matching coefficient signs. Opposite-sign significant matches are separately reported. Cohort-specific q values should not be substituted for the all-available genome-wide q values. Descriptive Spearman rho compares MLS partial-r profiles across genes. Nearby genes and overlapping windows are dependent, so no gene-independence-based significance claim is made for those profile correlations. There is no phylogenetic regression or independent experimental validation.

Figures: comparison_*_nominal.pdf/png shows rho plus positive/negative nominal overlaps. comparison_*_FDR.pdf/png shows the corresponding FDR overlaps. joint_correlation_map.pdf/png compares all seven methods on one common gene/species cohort. shared_gene_correlation_heatmap.pdf/png displays selected nominally shared genes; '*' means p<0.05, not FDR significance. The plot has no overall title so it can be composed with other panels.

Reproduction: use the project Python scientific environment (numpy, scipy, pyBigWig, matplotlib). Run scripts/quantify_activity.py, compare.py, plot.py and package_results.py from their original location under output/activity_omics_comparison_2026-10-20. Cached values are keyed by current input file size/mtime and locus hash. No source files, browser defaults or source measurements were modified.
'''
    (OUT/'METHODS.txt').write_text(methods)
    provenance={'activity':[json.loads(p.read_text()) for p in sorted((HERE/'activity_cache').glob('*.json'))]}
    for name in ('gene_activity.sqlite','liver_omics.sqlite','species_mls.json'):
        p=ROOT/'cross_species_visual_tool/web_app/data'/name
        provenance[name]={'path':str(p),'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest() if hasattr(hashlib,'file_digest') else hashlib.sha256(p.read_bytes()).hexdigest()}
    (OUT/'input_provenance.json').write_text(json.dumps(provenance,indent=2))
    (OUT/'scripts').mkdir(exist_ok=True)
    for p in HERE.glob('*.py'):shutil.copy2(p,OUT/'scripts'/p.name)
    archive=OUT.parent/'activity_omics_comparison_2026-10-20.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file():z.write(p,Path(OUT.name)/p.relative_to(OUT))
    print('\n'.join(report))
    print('Bundle:',archive)


if __name__=='__main__':main()
