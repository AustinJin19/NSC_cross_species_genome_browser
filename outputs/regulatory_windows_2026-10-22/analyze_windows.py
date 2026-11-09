"""Regulatory peak counts versus MLS at six TSS flanks.

Uses the browser's peak parser, chromosome resolver and 50%-span rule.
Vectorized counts are checked against the existing scorer and the saved 500 kb run.
"""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, t

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
APP = ROOT / 'cross_species_visual_tool/web_app'
sys.path.insert(0, str(APP))
from activity import DB
from server import SPECS, peak_path, resolve_chrom
from regulatory_strength import peak_arrays, score_records
from region_availability import chrom_sizes

WINDOWS = [20, 50, 100, 200, 500, 1000]
EXCLUDED = {'Human', 'Macaque', 'Balaena_mysticetus'}
LIFE = json.loads((DB.parent / 'species_mls.json').read_text())
FIELDS = ['window_kb', 'gene', 'n_species', 'partial_r', 'p_value',
          'MLS_coefficient', 'mass_coefficient_per_log10g', 'model_intercept',
          'df', 'species', 'BH_FDR_window', 'BH_FDR_all_windows']
(OUT / 'selected_loci').mkdir(parents=True, exist_ok=True)


def write(path, rows, fields):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def quantify(records, sid, folder):
    """Return (window, locus) medians and minimum measured-sample span fractions."""
    n = len(records)
    samples = sorted({s for r in records for s in r['samples']})
    counts = np.full((len(WINDOWS), len(samples), n), np.nan)
    fractions = np.full_like(counts, np.nan)
    chroms, starts = [], []
    for r in records:
        c, span, strand = r['locus'].rsplit('|', 1)[1].rsplit(':', 2)
        a, b = map(int, span.split('-'))
        chroms.append(c)
        starts.append(a if strand == '+' else b - 1)
    tss = np.asarray(starts, dtype=np.int64)
    grouped = defaultdict(list)
    for i, c in enumerate(chroms):
        grouped[c].append(i)
    for j, sample in enumerate(samples):
        bw = folder / (sample + '.bw')
        bed = peak_path(bw)
        if not bw.exists() or not bed.exists():
            continue
        stat = bw.stat()
        sizes = chrom_sizes(str(bw), stat.st_size, stat.st_mtime_ns)
        stat = bed.stat()
        peaks = peak_arrays(str(bed), stat.st_size, stat.st_mtime_ns)
        compatible = not peaks or any(resolve_chrom(c, sizes) for c in peaks)
        sample_present = np.array([sample in r['samples'] for r in records])
        for chrom, indexes in grouped.items():
            c = resolve_chrom(chrom, sizes)
            if c is None:
                continue
            pc = resolve_chrom(chrom, peaks)
            if pc is None and not compatible:
                continue
            ix = np.asarray(indexes, dtype=np.int64)
            ix = ix[sample_present[ix]]
            for k, kb in enumerate(WINDOWS):
                flank = kb * 1000
                lo = np.maximum(0, tss[ix] - flank)
                hi = np.minimum(sizes[c], tss[ix] + flank)
                width = np.maximum(0, hi - lo)
                valid = width > 0
                jj = ix[valid]
                value = (np.zeros(len(jj), dtype=np.int64) if pc is None else
                         np.searchsorted(peaks[pc][0], hi[valid], side='left') -
                         np.searchsorted(peaks[pc][1], lo[valid], side='right'))
                counts[k, j, jj] = value
                fractions[k, j, jj] = width[valid] / (2 * flank)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        medians = np.nanmedian(counts, axis=1)
        available = np.nanmin(fractions, axis=1)
    n_samples = np.isfinite(counts).sum(axis=1)
    # Cross-check deterministic loci, chromosome edges, duplicate loci and SIRT6
    # against the production browser scorer for every requested window.
    ix = set(np.linspace(0, max(n-1, 0), min(n, 10), dtype=int))
    ix.update(i for i, r in enumerate(records) if r['symbol_key'] in ('SIRT6', 'IGF1', 'HAS2'))
    ix = sorted(ix)
    for k, kb in enumerate(WINDOWS):
        reference = score_records([records[i] for i in ix], kb*1000)
        for i, ref in zip(ix, reference):
            observed = medians[k, i]
            if ref['median'] is None:
                assert np.isnan(observed), (sid, kb, ref['locus'])
            else:
                assert observed == ref['median'], (sid, kb, ref['locus'], observed, ref['median'])
                assert math.isclose(available[k, i], ref['available_fraction'], abs_tol=1e-15)
                assert bool(available[k, i] >= .5) == ref['region_eligible']
    return medians, available, n_samples, len(ix)*len(WINDOWS)


def fit(gene, values, traits, kb):
    ids = sorted(s for s, v in values.items() if v['region_eligible'] and v['median'] > 0)
    n = len(ids)
    if n < 5:
        return None
    x = np.array([math.log10(traits[s]['MLS_years']) for s in ids])
    y = np.array([math.log10(values[s]['median']) for s in ids])
    mass = np.array([math.log10(traits[s]['body_mass_g']) for s in ids])
    C = np.column_stack([np.ones(n), mass])
    if np.linalg.matrix_rank(C) < 2:
        return None
    rx = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
    ry = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
    if (rx @ rx <= 1e-12 * np.sum((x-x.mean())**2) or
            ry @ ry <= 1e-12 * np.sum((y-y.mean())**2)):
        return None
    r = float(np.clip(rx @ ry / math.sqrt((rx @ rx)*(ry @ ry)), -1, 1))
    df = n - 3
    p = float(2*t.sf(abs(r)*math.sqrt(df / max(1-r*r, 1e-300)), df))
    X = np.column_stack([np.ones(n), x, mass])
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y - X @ coef
    covariance = (residual @ residual)/df * np.linalg.inv(X.T @ X)
    p_ols = float(2*t.sf(abs(coef[1]/math.sqrt(covariance[1, 1])), df))
    assert abs(p - p_ols) < 1e-8, (gene, kb, p, p_ols)
    return dict(window_kb=kb, gene=gene, n_species=n, partial_r=r, p_value=p,
                MLS_coefficient=float(coef[1]), mass_coefficient_per_log10g=float(coef[2]),
                model_intercept=float(coef[0]), df=df, species=','.join(ids),
                ols_p_difference=abs(p-p_ols))


def main():
    started = time.time()
    genes = {kb: defaultdict(dict) for kb in WINDOWS}
    traits, source_files, coverage, species_rows = {}, [], [], []
    validation_count = 0
    con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    for sid, common, latin, color, folder, prefix in SPECS:
        trait = LIFE.get(latin, {})
        years, mass = trait.get('years'), trait.get('body_mass_g')
        reason = ('Excluded by analysis settings' if sid in EXCLUDED else
                  'MLS unavailable or lower bound' if not isinstance(years, (int, float)) or years <= 0 else
                  'Body mass unavailable' if not isinstance(mass, (int, float)) or mass <= 0 else 'Eligible')
        species_rows.append(dict(species=sid, common_name=common, scientific_name=latin,
                                 MLS_years=years, body_mass_g=mass, analysis_status=reason,
                                 trait_source=trait.get('source', ''), trait_note=trait.get('note', '')))
        if reason != 'Eligible':
            continue
        traits[sid] = species_rows[-1]
        records = [dict(r) for table in ['activity', 'activity_additions'] if table in tables
                   for r in con.execute(f'SELECT * FROM {table} WHERE species=? ORDER BY locus', (sid,))]
        for r in records:
            r['samples'] = json.loads(r['samples'])
        # Stable lexical locus order is the browser's tie-breaking convention.
        records.sort(key=lambda r: r['locus'])
        sample_folder = ROOT / 'cross_species_visual_tool' / folder
        medians, available, n_samples, n_check = quantify(records, sid, sample_folder)
        validation_count += n_check
        for sample in sorted({s for r in records for s in r['samples']}):
            bw = sample_folder/(sample+'.bw')
            for kind, path in [('bigWig', bw), ('peaks', peak_path(bw))]:
                st = path.stat() if path.exists() else None
                source_files.append(dict(species=sid, sample=sample, kind=kind, path=str(path),
                                         size=st.st_size if st else None, mtime_ns=st.st_mtime_ns if st else None))
        for k, kb in enumerate(WINDOWS):
            selected = {}
            for i, r in enumerate(records):
                median = medians[k, i]
                if not math.isfinite(median):
                    continue
                gene = r['symbol_key']
                if gene not in selected or median > selected[gene]['median']:
                    selected[gene] = dict(gene=gene, species=sid, locus=r['locus'],
                                          median=float(median), available_fraction=float(available[k, i]),
                                          region_eligible=bool(available[k, i] >= .5), n_samples=int(n_samples[k, i]))
            write(OUT/'selected_loci'/f'{sid}_{kb}kb.tsv.gz', selected.values(),
                  ['gene', 'species', 'locus', 'median', 'available_fraction', 'region_eligible', 'n_samples'])
            for gene, r in selected.items():
                genes[kb][gene][sid] = r
            coverage.append(dict(window_kb=kb, species=sid, loci=len(records), selected_genes=len(selected),
                                 excluded_region=sum(not r['region_eligible'] for r in selected.values()),
                                 excluded_nonpositive=sum(r['median'] <= 0 for r in selected.values())))
        print(f'Quantified all six windows: {sid}, {len(records):,} loci; elapsed {time.time()-started:.1f}s', flush=True)
    con.close()
    results = {}
    for kb in WINDOWS:
        rows = [result for gene, sp in genes[kb].items() if (result := fit(gene, sp, traits, kb)) is not None]
        q = false_discovery_control([r['p_value'] for r in rows], method='bh')
        for row, value in zip(rows, q):
            row['BH_FDR_window'] = float(value)
        results[kb] = sorted(rows, key=lambda r: (r['p_value'], r['gene']))
        print(f'Fitted ±{kb} kb: {len(rows):,} genes', flush=True)
    all_rows = [r for rows in results.values() for r in rows]
    q = false_discovery_control([r['p_value'] for r in all_rows], method='bh')
    for row, value in zip(all_rows, q):
        row['BH_FDR_all_windows'] = float(value)
    # The original 71-gene run must reproduce from unchanged source data.
    prior_path = ROOT/'outputs/regulatory-500kb-mass-2026-10-19/all_tested.tsv'
    prior = {r['gene']: r for r in csv.DictReader(prior_path.open(), delimiter='\t')}
    current = {r['gene']: r for r in results[500]}
    assert prior.keys() == current.keys(), '500 kb tested-gene set changed'
    maximum_difference = 0.0
    for gene, old in prior.items():
        new = current[gene]
        assert int(old['n_species']) == new['n_species'] and old['species'] == new['species']
        for a, b in [('partial_r', 'partial_r'), ('p_value', 'p_value'), ('BH_FDR', 'BH_FDR_window')]:
            difference = abs(float(old[a])-new[b])
            maximum_difference = max(maximum_difference, difference)
            assert difference < 1e-10, (gene, a, difference)
    summary = []
    significant = []
    for kb, rows in results.items():
        folder = OUT/f'tss_{kb}kb'
        folder.mkdir(exist_ok=True)
        write(folder/'all_tested.tsv.gz', rows, FIELDS)
        for name, column in [('FDR05', 'BH_FDR_window'), ('joint_FDR05', 'BH_FDR_all_windows'), ('nominal_p05', 'p_value')]:
            for direction, sign in [('positive', 1), ('negative', -1)]:
                sel = [r for r in rows if r[column] < .05 and r['partial_r']*sign > 0]
                write(folder/f'{name}_{direction}.tsv', sel, FIELDS)
                (folder/f'{name}_{direction}_genes.txt').write_text(''.join(r['gene']+'\n' for r in sel))
        sig = [dict(r, direction='Positive' if r['partial_r'] > 0 else 'Negative') for r in rows if r['BH_FDR_window'] < .05]
        significant.extend(sig)
        summary.append(dict(window_kb=kb, full_width_kb=2*kb, tested=len(rows),
                            FDR_positive=sum(r['partial_r'] > 0 for r in sig),
                            FDR_negative=sum(r['partial_r'] < 0 for r in sig),
                            joint_FDR_positive=sum(r['BH_FDR_all_windows'] < .05 and r['partial_r'] > 0 for r in rows),
                            joint_FDR_negative=sum(r['BH_FDR_all_windows'] < .05 and r['partial_r'] < 0 for r in rows),
                            nominal_positive=sum(r['p_value'] < .05 and r['partial_r'] > 0 for r in rows),
                            nominal_negative=sum(r['p_value'] < .05 and r['partial_r'] < 0 for r in rows),
                            minimum_FDR=min(r['BH_FDR_window'] for r in rows)))
    union = {r['gene'] for r in significant}
    matrix, wide = [], []
    lookup = {kb: {r['gene']: r for r in rows} for kb, rows in results.items()}
    for gene in union:
        sig_windows = [kb for kb in WINDOWS if gene in lookup[kb] and lookup[kb][gene]['BH_FDR_window'] < .05]
        record = dict(gene=gene, significant_windows=len(sig_windows), windows_kb='; '.join(map(str, sig_windows)))
        for kb in WINDOWS:
            r = lookup[kb].get(gene)
            record.update({f'r_{kb}kb': r['partial_r'] if r else None,
                           f'FDR_{kb}kb': r['BH_FDR_window'] if r else None,
                           f'n_{kb}kb': r['n_species'] if r else None})
        wide.append(record)
        for kb in WINDOWS:
            if gene in lookup[kb]:
                matrix.append(lookup[kb][gene])
    wide.sort(key=lambda r: (-r['significant_windows'], min(r[f'FDR_{kb}kb'] for kb in WINDOWS if r[f'FDR_{kb}kb'] is not None), r['gene']))
    overlap = []
    sets = {kb: {r['gene'] for r in significant if r['window_kb']==kb} for kb in WINDOWS}
    for kb in WINDOWS:
        overlap.append(dict(window_kb=kb, **{f'overlap_{w}kb':len(sets[kb]&sets[w]) for w in WINDOWS}))
    write(OUT/'window_summary.tsv', summary, list(summary[0]))
    write(OUT/'significant_genes_all_windows.tsv', significant, FIELDS+['direction'])
    write(OUT/'significant_gene_comparison.tsv', wide, list(wide[0]))
    write(OUT/'significant_gene_window_overlap.tsv', overlap, list(overlap[0]))
    write(OUT/'species.tsv', species_rows, list(species_rows[0]))
    write(OUT/'coverage_by_window_species.tsv', coverage, list(coverage[0]))
    report = dict(summary=summary, significant=significant, comparison=wide,
                  comparison_long=matrix, overlap=overlap, species=species_rows,
                  validation=dict(browser_scorer_locus_windows=validation_count,
                                  maximum_500kb_difference=maximum_difference,
                                  maximum_OLS_p_difference=max(r['ols_p_difference'] for r in all_rows)),
                  union_significant=len(union), all_six_significant=[r['gene'] for r in wide if r['significant_windows']==6],
                  sirt6=[lookup[kb].get('SIRT6') for kb in WINDOWS])
    (OUT/'report_data.json').write_text(json.dumps(report, indent=2))
    manifest = dict(windows_flank_kb=WINDOWS, measurement='median count of peaks overlapping TSS window',
                    model='log10(species median peak count) ~ log10(MLS years) + log10(body mass g)',
                    excluded_species=sorted(EXCLUDED), lower_bound_MLS_excluded=True,
                    minimum_species=5, minimum_available_fraction=.5,
                    locus_selection='highest species median before availability and positive-value filters; no fallback',
                    multiplicity='BH within each window; additional BH across all six gene-window tests',
                    availability='window overlap with bigWig chromosome/scaffold bounds; internal assembly gaps not measured',
                    species_eligible=len(traits), total_gene_window_tests=len(all_rows),
                    database=str(DB), database_mtime_ns=DB.stat().st_mtime_ns,
                    traits_sha256=hashlib.sha256((DB.parent/'species_mls.json').read_bytes()).hexdigest(),
                    source_files=source_files, validation=report['validation'])
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(dict(summary=summary, union_significant=len(union), all_six_significant=report['all_six_significant'], validation=report['validation']), indent=2), flush=True)


if __name__ == '__main__':
    main()
