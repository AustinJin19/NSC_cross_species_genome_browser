"""Verify selected results against current browser queries and statistics."""
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
import analyze as a

out = Path(__file__).resolve().parent
data = json.loads((out/'results.json').read_text())
rows = data['results']
traits = {r['species']: r for r in data['inventory'] if r['status']=='Included'}
selected = {'TERT', 'SIRT6', 'CREB3L3'}
for window in (2,5):
    for mass in ('Yes','No'):
        chosen = [r for r in rows if r['window_kb']==window and r['body_mass_correction']==mass]
        selected.update(r['gene'] for r in chosen[:2])
checks = []
details = defaultdict(dict)
with gzip.open(out/'selected_loci_and_coverage.tsv.gz','rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['gene'] in selected:
            details[(r['gene'],r['method'])][r['species']] = r
for gene in sorted(selected):
    for method,_ in a.METHODS:
        actual = a.activity.query(gene,method)
        grouped = defaultdict(list)
        for r in actual['records']:
            if r['species'] in traits and r['median'] is not None:
                grouped[r['species']].append(r)
        expected = details[(gene,method)]
        assert set(grouped)==set(expected), (gene,method,set(grouped)^set(expected))
        included = []
        max_delta = 0
        for sid, rr in grouped.items():
            r = sorted(rr,key=lambda r:(-r['median'],r['locus']))[0]
            e = expected[sid]
            assert r['locus']==e['locus'], (gene,method,sid,r['locus'],e['locus'])
            assert np.isclose(r['median'],float(e['median_activity']),rtol=2e-7,atol=1e-5), (gene,method,sid,r['median'],e['median_activity'])
            assert np.isclose(r['positive_signal_fraction'],float(e['positive_signal_fraction']),atol=1e-12), (gene,method,sid,'positive')
            assert np.isclose(r['available_fraction'],float(e['available_fraction']),atol=1e-12), (gene,method,sid,'available')
            max_delta=max(max_delta,abs(r['median']-float(e['median_activity'])))
            if r['median']>0 and r['positive_signal_fraction']>=.5 and r['available_fraction']>=.5:
                included.append(sid)
        for mass in ('Yes','No'):
            matches=[r for r in rows if r['gene']==gene and r['window_kb']==int(method[3:-2]) and r['body_mass_correction']==mass]
            if matches: assert sorted(matches[0]['species'].split('; '))==sorted(included)
        checks.append(dict(gene=gene,method=method,species_checked=len(grouped),species_in_fit=len(included),maximum_score_difference=max_delta))
        print(checks[-1],flush=True)
node_inputs=[]
for r in rows:
    if r['gene'] not in selected:continue
    obs=details[(r['gene'],'tss'+str(r['window_kb'])+'kb')]
    pts=[dict(x=float(np.log10(traits[sid]['MLS_years'])),y=float(np.log10(float(obs[sid]['median_activity']))),
              logMass=float(np.log10(traits[sid]['body_mass_g']))) for sid in r['species'].split('; ')]
    node_inputs.append(dict(gene=r['gene'],window=r['window_kb'],mass=r['body_mass_correction'],points=pts,r=r['r'],p=r['p'],slope=r['MLS_beta']))
(out/'browser_statistics_checks.json').write_text(json.dumps(node_inputs))
(out/'validation.json').write_text(json.dumps(dict(browser_query_checks=checks),indent=2))
print('Browser data and cohort verification passed.',flush=True)
