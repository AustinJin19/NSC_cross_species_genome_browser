#!/usr/bin/env python3
"""Compare MLS associations across ATAC windows, liver RNA and protein.

Primary comparisons use the same genes and, for each gene, exactly the same
species across all five ATAC methods and the selected expression assay.
"""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'cross_species_visual_tool/web_app').is_dir())
HERE = Path(__file__).resolve().parent
WEB = ROOT / 'cross_species_visual_tool/web_app'
sys.path[:0] = [str(WEB), str(WEB / '.deps'), str(HERE)]
import numpy as np
from scipy.stats import t, spearmanr, pearsonr, false_discovery_control
from quantify_activity import SELECTED, METHODS, CACHE
from omics import query

OUT = ROOT / 'outputs/activity_omics_comparison_2026-10-20'
OUT.mkdir(exist_ok=True, parents=True)
ATAC = [m for m, _ in METHODS]
ALL = ATAC + ['rna', 'protein']
LIFE = json.loads((WEB / 'data/species_mls.json').read_text())
SPECIES = {s[0]: s for s in SELECTED}
# This taxonomic name mapping already exists in the project's trait metadata.
ALIASES = {'Ellobius lutescens': 'Bramus lutescens'}
SID = {s[2]: s[0] for s in SELECTED}
for source, target in ALIASES.items():
    SID[source] = SID[target]
X = {sid: math.log10(LIFE[s[2]]['years']) for sid, s in SPECIES.items()}
MASS = {sid: math.log10(LIFE[s[2]]['body_mass_g']) for sid, s in SPECIES.items()}
LABELS = dict(tss2kb='TSS ±2 kb', tss5kb='TSS ±5 kb', body='Gene body + 2 kb',
              tss200kb='TSS ±200 kb', tss500kb='TSS ±500 kb', rna='Transcriptome', protein='Proteome')


def write(path, rows, fields=None):
    path.parent.mkdir(exist_ok=True, parents=True)
    if fields is None:
        fields = list(rows[0]) if rows else ['gene']
    with path.open('w') as h:
        writer = csv.DictWriter(h, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def load_inputs():
    data = {m: defaultdict(dict) for m in ALL}
    lengths = defaultdict(dict)
    excluded_region = defaultdict(int)
    for sid in SPECIES:
        with gzip.open(CACHE / f'{sid}.tsv.gz', 'rt') as h:
            for row in csv.DictReader(h, delimiter='\t'):
                m, g = row['method'], row['gene']
                if row['region_eligible'] != '1':
                    excluded_region[m] += 1
                    continue
                activity, length = float(row['median_activity']), float(row['length_kb'])
                if activity <= 0 or (m == 'body' and length <= 0):
                    continue
                data[m][g][sid] = math.log10(activity)
                if m == 'body':
                    lengths[g][sid] = math.log10(length)
    assert not query('SIRT6')['stale'], 'The RNA/protein index is stale'
    records = defaultdict(list)
    notes = set()
    with sqlite3.connect(f'file:{WEB}/data/liver_omics.sqlite?mode=ro', uri=True) as c:
        rows = c.execute('SELECT assay,symbol,species,sample,log,normalized,common,norm_on FROM measurements '
                         'WHERE NOT (assay="rna" AND sample="BMR_Turk_Lu")')
        seen = set()
        for assay, gene, species, sample, log, normalized, common, norm_on in rows:
            if species not in SID:
                continue
            key = assay, gene, SID[species], sample
            assert key not in seen, ('duplicate sample measurement', key)
            seen.add(key)
            val = log if assay == 'rna' else normalized
            if val is None or not math.isfinite(val):
                continue
            if assay == 'protein':
                assert norm_on == 'log.quant', norm_on
            records[(assay, gene, SID[species])].append(val)
            if species == 'Acomys cahirinus' and common == 'Acomys russatus':
                notes.add('RNA Acomys cahirinus species field has Acomys russatus in common field; ASM-excluded sensitivity supplied.')
        sources = json.loads(c.execute("SELECT value FROM metadata WHERE key='sources'").fetchone()[0])
    for (assay, gene, sid), values in records.items():
        data[assay][gene][sid] = float(np.median(values))
    write(OUT / 'species_inclusion.tsv', [dict(species=sid, common=s[1], scientific_name=s[2],
          MLS_years=LIFE[s[2]]['years'], body_mass_g=LIFE[s[2]]['body_mass_g'],
          RNA_species_present=int(any(sid in v for v in data['rna'].values())),
          protein_species_present=int(any(sid in v for v in data['protein'].values())),
          trait_source=LIFE[s[2]].get('source'), trait_note=LIFE[s[2]].get('note'))
          for sid, s in SPECIES.items()])
    with gzip.open(OUT / 'analysis_input_species_medians.tsv.gz', 'wt') as h:
        w = csv.writer(h, delimiter='\t')
        w.writerow(['method', 'gene', 'species', 'analysis_value', 'log10_gene_length_kb'])
        for method in ALL:
            for gene in sorted(data[method]):
                for sid, val in sorted(data[method][gene].items()):
                    w.writerow([method, gene, sid, val, lengths[gene].get(sid, '') if method == 'body' else ''])
    settings = dict(activity='log10(species median integrated CPM·bp); no pseudocount',
                    regions=dict(METHODS), MLS='log10(years)',
                    RNA='species median of source log.tpm; no further logarithm',
                    protein='species median of source quant.norm normalized on log.quant; no further logarithm',
                    primary_covariates='log10 body mass in all assays; additional log10 gene-body length for gene-body ATAC only',
                    species_exclusions=['Human', 'Rhesus macaque', 'Bowhead whale', 'missing or lower-bound MLS'],
                    sample_exclusion='Browser sample exclusions retained; BMR_Turk_Lu removed from RNA',
                    min_species=5, available_fraction_minimum=.5,
                    locus='Highest species-median locus separately for each method; choose before eligibility; no fallback',
                    primary_matching='For each RNA/protein comparison, intersection of valid species across all five ATAC methods and the expression assay for each gene; retain genes with valid fits for all six methods',
                    multiple_testing='BH separately for each method within each declared cohort and sensitivity analysis; nominal p<0.05 and BH q<0.05 reported separately',
                    joint_cohort='All seven methods use identical gene/species sets per gene',
                    matching='Uppercase gene symbol, not verified one-to-one orthology',
                    phylogenetic_correction=False, taxonomic_aliases=ALIASES,
                    notes=sorted(notes), omics_sources=sources, region_exclusions_species_gene=dict(excluded_region))
    (OUT / 'settings.json').write_text(json.dumps(settings, indent=2))
    return data, lengths, settings


def fit(x, y, mass=None, length=None):
    n = len(x)
    if n < 5 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    columns = [np.ones(n)]
    if mass is not None:
        columns.append(mass)
    if length is not None:
        columns.append(length)
    C = np.column_stack(columns)
    if np.linalg.matrix_rank(C) < C.shape[1]:
        return None
    Q, _ = np.linalg.qr(C, mode='reduced')
    rx, ry = x - Q @ (Q.T @ x), y - Q @ (Q.T @ y)
    xx, yy = rx @ rx, ry @ ry
    if xx <= 1e-12 * np.sum((x-x.mean())**2) or yy <= 1e-12 * np.sum((y-y.mean())**2):
        return None
    r = float(np.clip((rx @ ry)/math.sqrt(xx*yy), -1, 1))
    df = n-C.shape[1]-1
    if df < 1:
        return None
    p = float(2*t.sf(abs(r)*math.sqrt(df/max(1-r*r, 1e-300)), df))
    return dict(n_species=n, r=r, p=p, df=df, MLS_beta=float((rx @ ry)/xx))


def bh(rows):
    if rows:
        for row, q in zip(rows, false_discovery_control([r['p'] for r in rows])):
            row['q'] = float(q)


def fit_cohort(data, lengths, methods, mass_correct=True, exclude_asm=False, matched=True, min_species=5):
    genes = sorted(set.intersection(*(set(data[m]) for m in methods))) if matched else sorted(set.union(*(set(data[m]) for m in methods)))
    fits = {m: [] for m in methods}
    for gene in genes:
        common = set.intersection(*(set(data[m].get(gene, {})) for m in methods)) if matched else None
        gene_fits = {}
        for m in methods:
            ids = sorted(s for s in (common if matched else data[m].get(gene, {})) if not (exclude_asm and s == 'ASM'))
            if len(ids) < min_species:
                if matched:
                    break
                continue
            x = np.array([X[s] for s in ids]); y = np.array([data[m][gene][s] for s in ids])
            mass = np.array([MASS[s] for s in ids]) if mass_correct else None
            length = np.array([lengths[gene][s] for s in ids]) if m == 'body' else None
            result = fit(x, y, mass, length)
            if result is not None:
                gene_fits[m] = dict(gene=gene, **result, species=','.join(ids))
        if matched and len(gene_fits) != len(methods):
            continue
        for m, row in gene_fits.items():
            fits[m].append(row)
    for m in methods:
        bh(fits[m])
    return {m: {r['gene']: r for r in rows} for m, rows in fits.items()}


def pair_summary(a, b, label_a, label_b):
    genes = sorted(set(a) & set(b))
    ar, br = np.array([a[g]['r'] for g in genes]), np.array([b[g]['r'] for g in genes])
    if len(genes) < 3:
        return dict(method_a=label_a, method_b=label_b, tested_genes=len(genes))
    out = dict(method_a=label_a, method_b=label_b, tested_genes=len(genes),
               spearman_rho=float(spearmanr(ar, br).statistic), pearson_r=float(pearsonr(ar, br).statistic),
               same_direction=int(np.sum(ar*br > 0)), same_direction_percent=float(np.mean(ar*br > 0)*100),
               median_species=float(np.median([a[g]['n_species'] for g in genes])),
               minimum_species=min(a[g]['n_species'] for g in genes), maximum_species=max(a[g]['n_species'] for g in genes))
    for col, prefix in [('p', 'nominal'), ('q', 'FDR')]:
        sa, sb = {g for g in genes if a[g][col] < .05}, {g for g in genes if b[g][col] < .05}
        both = sa & sb
        pos = {g for g in both if a[g]['r'] > 0 and b[g]['r'] > 0}
        neg = {g for g in both if a[g]['r'] < 0 and b[g]['r'] < 0}
        out.update({f'{prefix}_a':len(sa), f'{prefix}_b':len(sb), f'{prefix}_both':len(both),
                    f'{prefix}_positive_both':len(pos), f'{prefix}_negative_both':len(neg),
                    f'{prefix}_same_direction':len(pos|neg), f'{prefix}_opposite_direction':len(both-pos-neg),
                    f'{prefix}_same_over_a_percent':len(pos|neg)/len(sa)*100 if sa else None,
                    f'{prefix}_jaccard':len(both)/len(sa|sb) if sa|sb else None})
    return out


def export_cohort(fits, variant, cohort):
    folder = OUT / variant / cohort
    combined = []
    for m, genes in fits.items():
        write(folder / f'{m}_correlations.tsv', sorted(genes.values(), key=lambda r:r['p']),
              ['gene', 'n_species', 'r', 'p', 'q', 'df', 'MLS_beta', 'species'])
        combined.extend(dict(method=m, **r) for r in genes.values())
    fields = ['method', 'gene', 'n_species', 'r', 'p', 'q', 'df', 'MLS_beta', 'species']
    write(folder / 'all_correlations.tsv', combined, fields)
    if cohort in ('rna', 'protein'):
        b = fits[cohort]
        summaries = []
        for m in ATAC:
            a = fits[m]
            summaries.append(dict(variant=variant, cohort=cohort, **pair_summary(a, b, m, cohort)))
            rows = []
            for gene in sorted(set(a) & set(b), key=lambda g:max(a[g]['p'], b[g]['p'])):
                ra, rb = a[gene], b[gene]
                direction = 'positive_both' if ra['r']>0 and rb['r']>0 else 'negative_both' if ra['r']<0 and rb['r']<0 else 'opposite'
                rows.append(dict(gene=gene, n_species=ra['n_species'], direction=direction,
                                 activity_r=ra['r'], activity_p=ra['p'], activity_q=ra['q'],
                                 expression_r=rb['r'], expression_p=rb['p'], expression_q=rb['q'], species=ra['species']))
            write(folder / f'{m}_vs_{cohort}_all.tsv', rows)
            for col, threshold in [('p','nominal_p05'), ('q','FDR05')]:
                for direction in ('positive_both', 'negative_both', 'opposite'):
                    selected = [r for r in rows if r[f'activity_{col}'] < .05 and r[f'expression_{col}'] < .05 and r['direction'] == direction]
                    write(folder / f'{m}_vs_{cohort}_{threshold}_{direction}.tsv', selected, list(rows[0]) if rows else None)
        return summaries
    return []


def verify_fits(fits, data, lengths, mass_correct):
    """Check partial-r p values using independent full-model OLS coefficient t tests."""
    validations = []
    for m, genes in fits.items():
        rows = sorted(genes.values(), key=lambda r:r['p'])
        indexes = sorted(set([0, len(rows)//2, len(rows)-1])) if rows else []
        for i in indexes:
            r = rows[i]; ids = r['species'].split(','); gene = r['gene']
            x = np.array([X[s] for s in ids]); y = np.array([data[m][gene][s] for s in ids])
            columns = [np.ones(len(ids)), x]
            if mass_correct: columns.append(np.array([MASS[s] for s in ids]))
            if m == 'body': columns.append(np.array([lengths[gene][s] for s in ids]))
            design = np.column_stack(columns)
            beta = np.linalg.lstsq(design, y, rcond=None)[0]
            rss = np.sum((y-design@beta)**2); df = len(y)-design.shape[1]
            se = math.sqrt((rss/df)*np.linalg.inv(design.T@design)[1,1])
            p = float(2*t.sf(abs(beta[1]/se), df))
            assert math.isclose(r['p'], p, rel_tol=1e-7, abs_tol=1e-10), (m, gene, p, r['p'])
            assert math.isclose(r['MLS_beta'], beta[1], rel_tol=1e-7, abs_tol=1e-10)
            validations.append(dict(method=m, gene=gene, p_partial_r=r['p'], p_full_OLS=p, df=df))
    return validations


def main():
    data, lengths, settings = load_inputs()
    summary, validation, method_counts = [], [], []
    for variant, mass, no_asm, minimum in [('mass_adjusted',True,False,5),
                                         ('without_mass',False,False,5),
                                         ('without_ASM',True,True,5),
                                         ('at_least_10_species',True,False,10)]:
        for cohort, methods in [('rna', ATAC+['rna']), ('protein', ATAC+['protein']), ('joint', ALL)]:
            fits = fit_cohort(data, lengths, methods, mass, no_asm, True, minimum)
            summary.extend(export_cohort(fits, variant, cohort))
            validation.extend(dict(variant=variant, cohort=cohort, **v) for v in verify_fits(fits, data, lengths, mass))
            for method, rows in fits.items():
                vals = list(rows.values())
                method_counts.append(dict(variant=variant, cohort=cohort, method=method, tested=len(vals),
                  nominal_positive=sum(r['p']<.05 and r['r']>0 for r in vals), nominal_negative=sum(r['p']<.05 and r['r']<0 for r in vals),
                  FDR_positive=sum(r['q']<.05 and r['r']>0 for r in vals), FDR_negative=sum(r['q']<.05 and r['r']<0 for r in vals)))
            print(variant, cohort, {m:len(v) for m,v in fits.items()}, flush=True)
            if variant == 'mass_adjusted' and cohort == 'joint':
                maps = [dict(cohort='joint', **pair_summary(fits[a], fits[b], a, b)) for a in ALL for b in ALL]
                write(OUT / 'joint_correlation_map.tsv', maps)
    # An all-available reference retains each method's own data coverage and
    # reproduces the latest genome-wide body analysis; it is NOT the matched heatmap.
    all_fits = fit_cohort(data, lengths, ALL, True, False, False)
    export_cohort(all_fits, 'mass_adjusted', 'all_available')
    validation.extend(dict(variant='mass_adjusted', cohort='all_available', **v) for v in verify_fits(all_fits, data, lengths, True))
    for method, rows in all_fits.items():
        vals=list(rows.values())
        method_counts.append(dict(variant='mass_adjusted',cohort='all_available',method=method,tested=len(vals),
          nominal_positive=sum(r['p']<.05 and r['r']>0 for r in vals),nominal_negative=sum(r['p']<.05 and r['r']<0 for r in vals),
          FDR_positive=sum(r['q']<.05 and r['r']>0 for r in vals),FDR_negative=sum(r['q']<.05 and r['r']<0 for r in vals)))
    write(OUT / 'comparison_summary.tsv', summary)
    write(OUT / 'method_significant_counts.tsv', method_counts)
    write(OUT / 'OLS_validation.tsv', validation)
    (OUT / 'summary.json').write_text(json.dumps(dict(comparisons=summary, counts=method_counts, settings=settings),indent=2))
    print(json.dumps([s for s in summary if s['variant']=='mass_adjusted'],indent=2),flush=True)


if __name__ == '__main__':
    main()
