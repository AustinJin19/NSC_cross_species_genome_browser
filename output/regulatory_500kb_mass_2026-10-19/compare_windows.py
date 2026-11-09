"""Validate genome-wide FDR and compare with the prior 50 kb analysis."""
import csv
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT = ROOT / 'outputs/regulatory-500kb-mass-2026-10-19'
new = json.loads((HERE / 'report_data.json').read_text())
old = json.loads((ROOT / 'output/regulatory_50kb_mass_2026-10-19/report_data.json').read_text())
rows = new['results']
p = np.array([r['p_value'] for r in rows])
assert np.all(p[:-1] <= p[1:])
q = np.minimum(1, np.minimum.accumulate((p * len(p) / np.arange(1, len(p)+1))[::-1])[::-1])
assert np.max(np.abs(q - np.array([r['BH_FDR'] for r in rows]))) < 1e-12
assert all(r['df'] == r['n_species']-3 for r in rows)
assert len({r['gene'] for r in rows}) == len(rows)
index = {r['gene']: r for r in rows}
comparison = []
for previous in old['results']:
    if previous['BH_FDR'] >= .05:
        continue
    current = index.get(previous['gene'])
    comparison.append(dict(gene=previous['gene'], partial_r_50kb=previous['partial_r'],
        FDR_50kb=previous['BH_FDR'], n_50kb=previous['n_species'],
        partial_r_500kb=current['partial_r'] if current else None,
        p_500kb=current['p_value'] if current else None,
        FDR_500kb=current['BH_FDR'] if current else None,
        n_500kb=current['n_species'] if current else None))
with (OUT/'comparison_previous_six.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(comparison[0]),delimiter='\t')
    w.writeheader();w.writerows(comparison)
(HERE/'comparison.json').write_text(json.dumps(comparison,indent=2))
print(json.dumps(comparison,indent=2))
for threshold,field in [('FDR05','BH_FDR'),('nominal_p05','p_value')]:
    for direction,sign in [('positive',1),('negative',-1)]:
        expected=[r['gene'] for r in rows if r[field]<.05 and r['partial_r']*sign>0]
        assert (OUT/f'{threshold}_{direction}_genes.txt').read_text().splitlines()==expected
print('Genome-wide BH, degrees of freedom and all gene-list exports verified.')
