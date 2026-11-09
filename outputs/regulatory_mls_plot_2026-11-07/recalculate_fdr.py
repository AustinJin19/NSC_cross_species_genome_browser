"""Current-browser +/-500 kb peak-count OLS tests for plot BH annotations.

Run with python3. Reuses the validated vectorized peak counter, but refreshes
samples, annotation/name mappings and traits. BH families are separate by model.
"""
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, t

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
APP = ROOT / 'cross_species_visual_tool/web_app'
sys.path[:0] = [str(APP), str(APP / '.deps')]
import activity
import gene_names
import server
from sample_qc import filter_records, QC_POLICY

spec = importlib.util.spec_from_file_location('window_counter', OUT.parent / 'regulatory_windows_2026-10-22/analyze_windows.py')
counter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(counter)
counter.WINDOWS = [500]
CACHE = OUT / 'fdr_quantification'
CACHE.mkdir(exist_ok=True)
EXCLUDED = {'Human', 'Macaque', 'Balaena_mysticetus'}
MIN_SPECIES = 5


def fingerprint(path):
    path = Path(path)
    st = path.stat() if path.exists() else None
    return dict(path=str(path.resolve()), size=st.st_size if st else None,
                mtime_ns=st.st_mtime_ns if st else None)


def write(path, rows):
    with gzip.open(path, 'wt', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def fit(gene, observations, traits, mass):
    ids = sorted(s for s, r in observations.items() if r['eligible'] and r['median'] > 0
                 and (not mass or isinstance(traits[s].get('body_mass_g'), (int,float))
                      and traits[s]['body_mass_g'] > 0))
    n = len(ids)
    if n < MIN_SPECIES: return None
    x = np.log10([traits[s]['years'] for s in ids])
    y = np.log10([observations[s]['median'] for s in ids])
    if np.ptp(x) == 0 or np.ptp(y) == 0: return None
    columns = [np.ones(n), x]
    if mass: columns.append(np.log10([traits[s]['body_mass_g'] for s in ids]))
    X = np.column_stack(columns)
    if np.linalg.matrix_rank(X) < X.shape[1]: return None
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    residual = y - X @ beta
    df = n - X.shape[1]
    se = math.sqrt((residual @ residual / df) * np.linalg.inv(X.T @ X)[1, 1])
    if not se > 0: return None
    tv = beta[1] / se
    p = float(2 * t.sf(abs(tv), df))
    return dict(gene=gene, strategy='peaks500kb', log10_activity=True, log10_MLS=True,
                body_mass_correction=mass, n_species=n, p_value=p,
                MLS_coefficient=float(beta[1]), standard_error=se,
                r=float(tv / math.sqrt(tv*tv + df)), df=df, species=';'.join(ids))


def main():
    started = time.time()
    life = json.loads((APP / 'data/species_mls.json').read_text())
    sources, traits, inputs = {}, {}, []
    resolver = gene_names.resolver()
    candidates = set()
    n_validated = 0
    with sqlite3.connect(f'file:{activity.DB}?mode=ro', uri=True) as con:
        con.row_factory = sqlite3.Row
        tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        for sid, common, latin, colour, folder, prefix in server.SPECS:
            trait = life.get(latin, {})
            if sid in EXCLUDED or not trait.get('years') or trait.get('lower_bound_years'): continue
            traits[sid] = trait
            raw = [dict(r) for table in ('activity','activity_additions') if table in tables
                   for r in con.execute(f'SELECT * FROM {table} WHERE species=?', (sid,))]
            records = filter_records(activity.refresh_sample_sources(raw))
            records = list({(r['symbol_key'],r['locus']):r for r in records}.values())
            records.sort(key=lambda r:r['locus'])
            folder = server.ROOT / folder
            samples = sorted({s for r in records for s in r['samples']})
            files = [p for sample in samples for p in (folder/(sample+'.bw'), server.peak_path(folder/(sample+'.bw')))]
            provenance = dict(species=sid, loci=[(r['symbol_key'],r['locus']) for r in records],
                              files=[fingerprint(p) for p in files], windows=[500], qc=QC_POLICY)
            digest = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
            cache = CACHE / (sid+'.npz')
            if cache.exists():
                with np.load(cache) as z:
                    if str(z['fingerprint']) != digest: raise ValueError('Stale quantification cache: '+sid)
                    medians, available = z['medians'], z['available']
            else:
                medians, available, n_samples, checked = counter.quantify(records, sid, folder)
                n_validated += checked
                np.savez_compressed(cache, medians=medians, available=available, fingerprint=digest)
            (CACHE/(sid+'.json')).write_text(json.dumps(provenance, indent=2))
            by_symbol = defaultdict(list)
            for i,r in enumerate(records):
                by_symbol[r['symbol_key']].append(i)
                candidates.add(resolver.canonical(r['symbol_key']))
            sources[sid] = (records, medians[0], available[0], by_symbol)
            inputs.extend(provenance['files'])
            print(f'{sid}: {len(records):,} loci; elapsed {time.time()-started:.1f}s', flush=True)
    families = {False:[], True:[]}
    for gene in sorted(candidates):
        observations = {}
        for sid,(records,medians,available,by_symbol) in sources.items():
            match = resolver.resolve(gene,sid,by_symbol)
            if match['symbol'] is None: continue
            indexes = [i for i in by_symbol[match['symbol']] if np.isfinite(medians[i])]
            if not indexes: continue
            i = min(indexes,key=lambda j:(-medians[j],records[j]['locus']))
            observations[sid] = dict(median=float(medians[i]),
                eligible=bool(np.isfinite(available[i]) and available[i]>=.5))
        if len(observations) < MIN_SPECIES: continue
        for mass in families:
            result = fit(gene,observations,traits,mass)
            if result is not None: families[mass].append(result)
    summary = []
    for mass, rows in families.items():
        q = false_discovery_control([r['p_value'] for r in rows], method='bh')
        for row,value in zip(rows,q): row['BH_FDR'] = float(value)
        rows.sort(key=lambda r:(r['p_value'],r['gene']))
        name = 'with_body_mass' if mass else 'without_body_mass'
        write(OUT/('peaks500kb_'+name+'_all_tested.tsv.gz'),rows)
        sirt6 = next(r for r in rows if r['gene']=='SIRT6')
        summary.append(dict(body_mass_correction=mass,genes_tested=len(rows),sirt6=sirt6))
        print(json.dumps(summary[-1]),flush=True)
    manifest = dict(window_kb=500,measurement='peak count',model='log10(median peak count) ~ log10(MLS) + optional log10(body mass)',
        excluded_species=sorted(EXCLUDED),macaque_included='Macaque' not in EXCLUDED,minimum_species=MIN_SPECIES,
        minimum_available_fraction=.5,positive_signal_filter=False,lower_bound_MLS_excluded=True,
        multiplicity='BH separately across all valid gene tests in each model; both signs included',
        naming='Current browser canonical names and resolver, ambiguous mappings omitted',
        locus_selection='Highest species median before filtering; lexical tie break',
        qc_policy=QC_POLICY,summary=summary,source_files=inputs,
        metadata_files=[fingerprint(activity.DB),fingerprint(gene_names.DATA),fingerprint(APP/'data/species_mls.json')],
        validation=dict(browser_scorer_locus_windows=n_validated),elapsed_seconds=time.time()-started)
    (OUT/'fdr_manifest.json').write_text(json.dumps(manifest,indent=2))


if __name__ == '__main__': main()
