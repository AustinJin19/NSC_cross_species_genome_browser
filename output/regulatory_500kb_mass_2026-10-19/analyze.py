"""Genome-wide regulatory peak count analysis matching current browser controls."""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, t

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
APP = ROOT / 'cross_species_visual_tool/web_app'
sys.path.insert(0, str(APP))
from activity import DB
from regulatory_strength import score_records
from server import SPECS, peak_path

OUT = ROOT / 'outputs/regulatory-500kb-mass-2026-10-19'
OUT.mkdir(parents=True, exist_ok=True)
(HERE / 'species_scores').mkdir(exist_ok=True)
LIFE = json.loads((DB.parent / 'species_mls.json').read_text())
FLANK = 500000
EXCLUDED = {'Human', 'Macaque', 'Balaena_mysticetus'}
selected = [s for s in SPECS if s[0] not in EXCLUDED
            and isinstance(LIFE.get(s[2], {}).get('years'), (int, float))
            and LIFE[s[2]]['years'] > 0 and LIFE[s[2]].get('body_mass_g', 0) > 0]
latin = {s[0]: s[2] for s in selected}
genes = defaultdict(dict)
provenance = []
con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
con.row_factory = sqlite3.Row
tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}


def write(path, rows, fields):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


for sid, name, scientific, color, folder, prefix in selected:
    started = time.time()
    records = [dict(r) for table in ['activity', 'activity_additions'] if table in tables
               for r in con.execute(f'SELECT * FROM {table} WHERE species=? ORDER BY locus', (sid,))]
    for r in records:
        r['samples'] = json.loads(r['samples'])
    scored = score_records(records, FLANK)
    # Preserve the same highest-median selection before any availability/log filters.
    for r in sorted(scored, key=lambda r: r['locus']):
        if r['median'] is None or not math.isfinite(r['median']):
            continue
        old = genes[r['symbol_key']].get(sid)
        if old is None or r['median'] > old['median']:
            genes[r['symbol_key']][sid] = r
    with gzip.open(HERE / 'species_scores' / f'{sid}.json.gz', 'wt') as f:
        json.dump(scored, f)
    source_files = []
    for sample in sorted({sample for r in records for sample in r['samples']}):
        bw = ROOT / 'cross_species_visual_tool' / folder / (sample + '.bw')
        for kind, path in [('bigWig', bw), ('peaks', peak_path(bw))]:
            stat = path.stat() if path.exists() else None
            source_files.append(dict(sample=sample, kind=kind, path=str(path),
                                     size=stat.st_size if stat else None,
                                     mtime_ns=stat.st_mtime_ns if stat else None))
    sp_selected = [sp[sid] for sp in genes.values() if sid in sp]
    info = dict(species=sid, scientific_name=scientific, MLS=LIFE[scientific]['years'],
                body_mass_g=LIFE[scientific]['body_mass_g'], loci=len(records),
                scored_loci=sum(r['median'] is not None for r in scored),
                selected_genes=len(sp_selected),
                excluded_region=sum(r['region_eligible'] is not True for r in sp_selected),
                excluded_zero=sum(r['median'] <= 0 for r in sp_selected), source_files=source_files)
    provenance.append(info)
    print(sid, 'loci', len(records), 'elapsed_s', round(time.time()-started, 1), flush=True)
con.close()
results = []
validation = {}
ols_error = 0.0
for gene, sp in genes.items():
    ids = sorted(sid for sid in sp if sp[sid]['region_eligible'] is True and sp[sid]['median'] > 0)
    n = len(ids)
    if n < 5:
        continue
    x = np.array([math.log10(LIFE[latin[sid]]['years']) for sid in ids])
    y = np.array([math.log10(sp[sid]['median']) for sid in ids])
    mass = np.array([math.log10(LIFE[latin[sid]]['body_mass_g']) for sid in ids])
    C = np.column_stack([np.ones(n), mass])
    if np.linalg.matrix_rank(C) < 2:
        continue
    rx = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
    ry = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
    if rx @ rx <= 1e-12 * np.sum((x-x.mean())**2) or ry @ ry <= 1e-12 * np.sum((y-y.mean())**2):
        continue
    r = float(np.clip((rx @ ry) / math.sqrt((rx @ rx)*(ry @ ry)), -1, 1))
    df = n - 3
    p = float(2*t.sf(abs(r)*math.sqrt(df / max(1-r*r, 1e-300)), df))
    X = np.column_stack([np.ones(n), x, mass])
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    # Independently compare the full-model slope t-test against the residual correlation.
    residual = y-X@coef
    covariance = (residual@residual)/df * np.linalg.inv(X.T@X)
    p_ols = float(2*t.sf(abs(coef[1]/math.sqrt(covariance[1, 1])), df))
    ols_error = max(ols_error, abs(p-p_ols))
    assert abs(p-p_ols) < 1e-8
    result = dict(gene=gene, n_species=n, partial_r=r, p_value=p,
                  MLS_coefficient=float(coef[1]), mass_coefficient_per_log10g=float(coef[2]),
                  model_intercept=float(coef[0]), df=df, species=','.join(ids))
    results.append(result)
    validation[gene] = dict(result=result, points=[dict(species=sid, locus=sp[sid]['locus'],
                            median=sp[sid]['median'], available_fraction=sp[sid]['available_fraction'],
                            x=float(x[i]), y=float(y[i]), logMass=float(mass[i])) for i, sid in enumerate(ids)])
q = false_discovery_control([r['p_value'] for r in results], method='bh')
for result, value in zip(results, q):
    result['BH_FDR'] = float(value)
results.sort(key=lambda r: r['p_value'])
fields = list(results[0])
write(OUT / 'all_tested.tsv', results, fields)
for threshold, column in [('FDR05', 'BH_FDR'), ('nominal_p05', 'p_value')]:
    for direction, sign in [('positive', 1), ('negative', -1)]:
        rows = [r for r in results if r[column] < .05 and r['partial_r']*sign > 0]
        write(OUT / f'{threshold}_{direction}.tsv', rows, fields)
        (OUT / f'{threshold}_{direction}_genes.txt').write_text(''.join(r['gene']+'\n' for r in rows))
summary = dict(tested=len(results), species_available=len(selected), minimum_species=5,
               FDR_positive=sum(r['BH_FDR'] < .05 and r['partial_r'] > 0 for r in results),
               FDR_negative=sum(r['BH_FDR'] < .05 and r['partial_r'] < 0 for r in results),
               nominal_positive=sum(r['p_value'] < .05 and r['partial_r'] > 0 for r in results),
               nominal_negative=sum(r['p_value'] < .05 and r['partial_r'] < 0 for r in results),
               minimum_FDR=min(r['BH_FDR'] for r in results), maximum_OLS_p_difference=ols_error,
               top_positive=[r for r in results if r['partial_r'] > 0][:10],
               top_negative=[r for r in results if r['partial_r'] < 0][:10])
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2))
(HERE / 'species_provenance.json').write_text(json.dumps(provenance, indent=2))
(HERE / 'report_data.json').write_text(json.dumps(dict(summary=summary, results=results, species=provenance)))
check_genes = list(dict.fromkeys(['USP9X', 'SIRT6', 'LMNA', results[0]['gene'], summary['top_positive'][0]['gene'], summary['top_negative'][0]['gene']]))
(HERE / 'validation_points.json').write_text(json.dumps({g: validation[g] for g in check_genes if g in validation}, indent=2))
settings = dict(measurement='peak count', flank_bp=FLANK, full_window_bp=2*FLANK,
                activity_scale='log10(species median peak count)', lifespan='log10(MLS)',
                covariates=['log10(body mass in grams)'], gene_length_correction=False,
                species_scope='all', exclude_species=sorted(EXCLUDED), min_available_fraction=.5,
                lower_bound_MLS_excluded=True, minimum_species=5, minimum_positive_activity=0,
                selection='highest median locus before availability filter; no fallback',
                database=str(DB), database_mtime_ns=DB.stat().st_mtime_ns,
                traits_sha256=hashlib.sha256((DB.parent/'species_mls.json').read_bytes()).hexdigest())
(HERE / 'settings.json').write_text(json.dumps(settings, indent=2))
print(json.dumps(summary, indent=2), flush=True)
