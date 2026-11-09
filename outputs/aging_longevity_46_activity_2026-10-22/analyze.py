"""Add validated gene-activity MLS associations to the 46 regulatory overlaps."""
import csv
import gzip
import hashlib
import json
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'outputs/aging_longevity_reference_2026-10-22'
ANALYSIS = ROOT / 'outputs/activity_omics_comparison_2026-10-20'
sys.path.insert(0, str(ROOT / 'output/activity_omics_comparison_2026-10-20'))
import quantify_activity as qa
import numpy as np
from scipy.stats import t, false_discovery_control


def read(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter='\t'))


def write(name, rows):
    with (OUT / name).open('w') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


reference = read(SOURCE / 'regulatory_overlap_unique_genes.tsv')
genes = {r['human_gene'] for r in reference}
assert len(genes) == len(reference) == 46
methods = {'tss2kb': 'TSS ±2 kb', 'body': 'Gene body + 2 kb promoter'}
specs = {s[0]: s for s in qa.SELECTED}
audit = []
cache = {}
for sid, spec in specs.items():
    loci = qa.read_loci(sid)
    files = sorted(p for p in (ROOT / 'cross_species_visual_tool' / spec[4]).glob('*.bw')
                   if p.name.lower().startswith(spec[5].lower()) and
                   p.stem.removesuffix('.cpm').removesuffix('_scaled') not in qa.EXCLUDED_SAMPLES)
    fingerprint_inputs = dict(version=1, loci_sha256=hashlib.sha256(json.dumps(loci).encode()).hexdigest(),
                              methods=qa.METHODS,
                              files=[dict(path=str(p.resolve()), size=p.stat().st_size,
                                          mtime_ns=p.stat().st_mtime_ns) for p in files])
    fingerprint = hashlib.sha256(json.dumps(fingerprint_inputs, sort_keys=True).encode()).hexdigest()
    meta = json.loads((qa.CACHE / f'{sid}.json').read_text())
    assert fingerprint == meta['fingerprint'], ('Stale activity cache', sid)
    audit.append(dict(species=sid, cache_fingerprint=fingerprint, matched_current_inputs=True))
    with gzip.open(qa.CACHE / f'{sid}.tsv.gz', 'rt') as h:
        for row in csv.DictReader(h, delimiter='\t'):
            if row['gene'] in genes and row['method'] in methods:
                cache[row['method'], row['gene'], sid] = row
print(f'Confirmed current annotation and bigWig fingerprints for {len(specs)} species.', flush=True)

old_traits = {r['species']: r for r in read(ANALYSIS / 'species_inclusion.tsv')}
assert set(old_traits) == set(specs)
for sid, s in specs.items():
    traits = qa.LIFE[s[2]]
    assert float(old_traits[sid]['MLS_years']) == traits['years'], sid
    assert float(old_traits[sid]['body_mass_g']) == traits['body_mass_g'], sid

tables, checks, long, species_rows = {}, [], [], []
for method in methods:
    rows = read(ANALYSIS / f'mass_adjusted/all_available/{method}_correlations.tsv')
    q = false_discovery_control([float(r['p']) for r in rows])
    assert np.allclose(q, [float(r['q']) for r in rows], rtol=1e-12, atol=1e-14)
    tables[method] = {r['gene']: r for r in rows}
    assert genes <= set(tables[method])
    for gene in sorted(genes):
        saved = tables[method][gene]
        included = []
        for sid, spec in specs.items():
            row = cache.get((method, gene, sid))
            trait = qa.LIFE[spec[2]]
            a = float(row['median_activity']) if row else None
            length = float(row['length_kb']) if row else None
            status = ('No activity measurement' if row is None else
                      'Available region below 50%' if row['region_eligible'] != '1' else
                      'Nonpositive activity' if a <= 0 else
                      'Invalid gene length' if method == 'body' and length <= 0 else 'Included')
            species_rows.append(dict(Gene=gene, Method=methods[method], Species=spec[1],
                Species_ID=sid, Scientific_name=spec[2], Median_activity_CPM_bp=a,
                Log10_median_activity=math.log10(a) if a is not None and a > 0 else None,
                MLS_years=trait['years'], Body_mass_g=trait['body_mass_g'], Gene_length_kb=length,
                Available_fraction=float(row['available_fraction']) if row else None,
                Sample_count=int(row['n_samples']) if row else None, Status=status,
                Locus=row['locus'] if row else None))
            if status == 'Included':
                included.append(sid)
        ids = sorted(included)
        assert ids == saved['species'].split(','), (gene, method, 'species changed')
        x = np.array([math.log10(qa.LIFE[specs[s][2]]['years']) for s in ids])
        y = np.array([math.log10(float(cache[method, gene, s]['median_activity'])) for s in ids])
        mass = np.array([math.log10(qa.LIFE[specs[s][2]]['body_mass_g']) for s in ids])
        columns = [np.ones(len(ids)), x, mass]
        if method == 'body':
            columns.append(np.array([math.log10(float(cache[method, gene, s]['length_kb'])) for s in ids]))
        design = np.column_stack(columns)
        beta, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
        assert rank == design.shape[1]
        df = len(y) - int(rank)
        rss = np.sum((y - design @ beta)**2)
        se = math.sqrt((rss / df) * np.linalg.inv(design.T @ design)[1, 1])
        statistic = beta[1] / se
        p = float(2 * t.sf(abs(statistic), df))
        r = float(statistic / math.sqrt(statistic**2 + df))
        for name, value in [('p', p), ('r', r), ('MLS_beta', float(beta[1]))]:
            assert math.isclose(value, float(saved[name]), rel_tol=1e-7, abs_tol=1e-10), (method, gene, name)
        checks.append(dict(method=method, gene=gene, p_from_OLS=p, p_from_source=float(saved['p']),
                           r_from_OLS=r, n_species=len(ids), df=df))
        long.append(dict(Gene=gene, Method=methods[method], Direction='Positive' if r > 0 else 'Negative',
            Partial_r=r, P=p, Genomewide_BH_FDR=float(saved['q']), N_species=len(ids),
            Nominal_P_below_005=p < .05, FDR_below_005=float(saved['q']) < .05,
            MLS_beta=float(beta[1]), DF=df, Species_IDs=','.join(ids)))

by_key = {(r['Gene'], r['Method']): r for r in long}
wide = []
for ref in reference:
    g = ref['human_gene']
    row = dict(Gene=g, Resources=ref['Data resources'],
               Regulatory_significant_windows_kb=ref['significant_windows_kb'],
               Regulatory_MLS_direction=ref['MLS_direction'])
    for m, label in methods.items():
        result = by_key[g, label]
        for col in ['Direction', 'Partial_r', 'P', 'Genomewide_BH_FDR', 'N_species',
                    'Nominal_P_below_005', 'FDR_below_005', 'MLS_beta']:
            row[f'{m}_{col}'] = result[col]
    row['Activity_direction_agreement'] = ('Same' if row['tss2kb_Direction'] == row['body_Direction'] else 'Opposite')
    row['Gene_name'] = ref['Gene name']
    for k in ['Experiment', 'GenAge', 'CellAge', 'GWAS', 'GenAge human', 'GenAge animal homolog',
              'GenAge fungal homolog', 'CellAge manipulation', 'CellAge expression',
              'GWAS longevity', 'GWAS healthspan', 'GWAS biological age', 'GWAS combined traits']:
        row[k] = ref[k]
    wide.append(row)

summary = []
for m, label in methods.items():
    selected = [r for r in long if r['Method'] == label]
    summary.append(dict(Method=label, Genes=len(selected),
        Positive=sum(r['Partial_r'] > 0 for r in selected), Negative=sum(r['Partial_r'] < 0 for r in selected),
        Nominal_positive=sum(r['P'] < .05 and r['Partial_r'] > 0 for r in selected),
        Nominal_negative=sum(r['P'] < .05 and r['Partial_r'] < 0 for r in selected),
        FDR_positive=sum(r['Genomewide_BH_FDR'] < .05 and r['Partial_r'] > 0 for r in selected),
        FDR_negative=sum(r['Genomewide_BH_FDR'] < .05 and r['Partial_r'] < 0 for r in selected),
        Genomewide_tests=len(tables[m]), Min_species=min(r['N_species'] for r in selected),
        Max_species=max(r['N_species'] for r in selected)))
direction_counts = Counter(f"{r['tss2kb_Direction']} / {r['body_Direction']}" for r in wide)
same_species = sum(tables['tss2kb'][g]['species'] == tables['body'][g]['species'] for g in genes)
write('46_genes_with_activity_MLS.tsv', wide)
write('correlations_long.tsv', long)
write('species_activity_values.tsv', species_rows)
write('summary.tsv', summary)

notes = [
    dict(Item='Scope', Details='The same 46 aging/longevity-reference overlaps from the six-window regulatory peak-count search. All 46 had positive regulatory-activity MLS associations in at least one window.'),
    dict(Item='TSS ±2 kb', Details='Unweighted integrated bigWig CPM signal over 2 kb upstream and 2 kb downstream of the annotated gene-locus TSS (4 kb total), clipped at contig boundaries. Not restricted to peaks.'),
    dict(Item='Gene body + 2 kb promoter', Details='Unweighted integrated bigWig CPM signal over the complete gene body plus 2 kb upstream in the strand-aware direction. Not restricted to peaks.'),
    dict(Item='Aggregation', Details='Median across available samples within species; choose the locus with the highest species median separately for each method before applying eligibility; no fallback to another locus.'),
    dict(Item='Model', Details='Response = log10(species median CPM·bp), predictor = log10(MLS in years). Both models adjust for log10(body mass in g). Gene-body model additionally adjusts for log10(gene-body length in kb). Partial r and two-sided OLS t-test P are reported.'),
    dict(Item='Exclusions', Details='Human, rhesus macaque and bowhead whale excluded. Missing/non-exact MLS or missing body mass, nonpositive activity and <50% available genomic region excluded. At least 5 species per fit. Available fraction concerns contig boundaries, not ATAC coverage.'),
    dict(Item='Species sets', Details=f'Each method uses all eligible species for that gene, without requiring RNA/protein data. {same_species}/46 genes have identical species sets between the two methods; see per-gene N and Species values for differences.'),
    dict(Item='Significance', Details='P < 0.05 is nominal only. BH FDR is retained from all 20,244 tested genes for TSS ±2 kb and all 20,200 for gene body; not recalculated only over the selected 46. No gene in these 46 has FDR < 0.05 for either activity method.'),
    dict(Item='Interpretation', Details='Positive/negative describes a cross-species association, not an experimentally established effect of gene activation or inhibition. Activity is integrated accessibility, not RNA expression. Symbol matching is not verified one-to-one orthology; no phylogenetic correction. Genes were selected using related ATAC data, so this is exploratory follow-up rather than independent replication.'),
    dict(Item='Verification', Details=f'Validated current annotations and bigWig size/mtime fingerprints for {len(specs)} species, current trait values, all 92 fits against independent full-model OLS, and genome-wide BH values.'),
    dict(Item='46-gene source', Details=str(SOURCE / 'regulatory_overlap_unique_genes.tsv')),
    dict(Item='Correlation source', Details=str(ANALYSIS / 'mass_adjusted/all_available')),
    dict(Item='Species activity source', Details=str(qa.CACHE)),
    dict(Item='Trait source', Details=str(ROOT / 'cross_species_visual_tool/web_app/data/species_mls.json')),
]
write('methods.tsv', notes)
data = dict(summary=summary, wide=wide, long=long, species_values=species_rows, notes=notes,
            direction_counts=dict(direction_counts), same_species_genes=same_species,
            nominal=[r for r in long if r['P'] < .05])
(OUT / 'results.json').write_text(json.dumps(data, indent=2))
source_paths = [SOURCE / 'regulatory_overlap_unique_genes.tsv',
                ROOT / 'cross_species_visual_tool/web_app/data/species_mls.json'] + [
                ANALYSIS / f'mass_adjusted/all_available/{m}_correlations.tsv' for m in methods]
(OUT / 'validation.json').write_text(json.dumps(dict(input_checks=audit, OLS_checks=checks,
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}), indent=2))

os.environ.setdefault('MPLCONFIGDIR', '/tmp/aging46_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
plt.rcParams.update({'font.family':'Arial','pdf.fonttype':42,'ps.fonttype':42})
order = sorted(wide, key=lambda r: (r['tss2kb_Direction'] != 'Positive', r['body_Direction'] != 'Positive',
                                    -r['body_Partial_r'], r['Gene']))
fig, axes = plt.subplots(1, 2, figsize=(9.5, 8.3), gridspec_kw={'wspace':.65})
norm = TwoSlopeNorm(vmin=-.65, vcenter=0, vmax=.65)
for panel, ax in enumerate(axes):
    subset = order[panel*23:(panel+1)*23]
    matrix = np.array([[r['tss2kb_Partial_r'], r['body_Partial_r']] for r in subset])
    im = ax.imshow(matrix, aspect='auto', cmap='RdBu_r', norm=norm)
    ax.set_yticks(range(len(subset)), [r['Gene'] for r in subset], fontsize=10)
    ax.set_xticks([0, 1], ['TSS ±2 kb', 'Gene body\n+ 2 kb promoter'], fontsize=10)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', which='both', length=0, pad=7)
    for i, r in enumerate(subset):
        for j, m in enumerate(methods):
            value = r[f'{m}_Partial_r']
            mark = '*' if r[f'{m}_P'] < .05 else ''
            ax.text(j, i, f'{value:+.2f}{mark}', ha='center', va='center', fontsize=9,
                    color='white' if abs(value) > .43 else '#17232E')
    ax.set_xticks(np.arange(-.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(subset), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.7)
    for spine in ax.spines.values():
        spine.set_visible(False)
fig.subplots_adjust(left=.14, right=.96, top=.88, bottom=.18)
cax = fig.add_axes([.35, .115, .37, .018])
cb = fig.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('Partial r with log10(MLS)', fontsize=10)
cb.ax.tick_params(labelsize=9)
fig.suptitle('Gene activity–lifespan associations in the 46 overlapping genes', fontsize=13, y=.985)
fig.text(.5, .94, 'log10(median activity); body-mass adjusted; gene-body method also length adjusted', ha='center', fontsize=9)
fig.text(.5, .025, '* Nominal P < 0.05; none pass genome-wide BH FDR < 0.05.\nHuman, macaque and bowhead whale excluded; ≥50% available region.', ha='center', fontsize=9)
fig.savefig(OUT / '46_genes_activity_MLS_heatmap.png', dpi=220, facecolor='white')
fig.savefig(OUT / '46_genes_activity_MLS_heatmap.pdf', facecolor='white')
print(json.dumps(dict(summary=summary, direction_counts=dict(direction_counts),
                     same_species_genes=same_species, nominal=data['nominal']), indent=2))
