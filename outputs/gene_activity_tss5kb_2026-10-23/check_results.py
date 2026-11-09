"""Verify current TSS ±5 kb activity inputs and fit MLS with/without mass."""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, t

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
ACT=ROOT/'outputs/activity_omics_comparison_2026-10-20'
WEB=ROOT/'cross_species_visual_tool/web_app'
sys.path[:0]=[str(WEB),str(WEB/'.deps')]
from server import SPECS, EXCLUDED_SAMPLES

def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def write(path,rows,fields):
    with path.open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)

traits={r['species']:r for r in read(ACT/'species_inclusion.tsv')}
life=json.loads((WEB/'data/species_mls.json').read_text())
specs={s[0]:s for s in SPECS if s[0] not in ('Human','Macaque','Balaena_mysticetus')
       and life.get(s[2],{}).get('years') and life[s[2]].get('body_mass_g')}
assert set(traits)==set(specs)
provenance=json.loads((ACT/'input_provenance.json').read_text())['activity']
checks=[]
with sqlite3.connect(f'file:{WEB}/data/gene_activity.sqlite?mode=ro',uri=True) as db:
    for previous in provenance:
        sid=previous['species'];s=specs[sid]
        assert float(traits[sid]['MLS_years'])==life[s[2]]['years']
        assert float(traits[sid]['body_mass_g'])==life[s[2]]['body_mass_g']
        loci=[]
        for table in ('activity','activity_additions'):
            loci.extend(db.execute(f'SELECT symbol_key,locus FROM {table} WHERE species=?',(sid,)))
        loci=sorted(set(loci),key=lambda r:r[1])
        assert hashlib.sha256(json.dumps(loci).encode()).hexdigest()==previous['inputs']['loci_sha256'],sid
        files=sorted(p for p in (ROOT/'cross_species_visual_tool'/s[4]).glob('*.bw')
                     if p.name.lower().startswith(s[5].lower()) and
                     p.stem.removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES)
        actual=[dict(path=str(p.resolve()),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in files]
        assert actual==previous['inputs']['files'],sid
        checks.append(dict(species=sid,samples=len(files),loci=len(loci)))
print('Verified current activity sources and traits for',len(checks),'species.',flush=True)

data=defaultdict(dict)
with gzip.open(ACT/'analysis_input_species_medians.tsv.gz','rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['method']=='tss5kb':data[r['gene']][r['species']]=float(r['analysis_value'])
prior={r['gene']:r for r in read(ACT/'mass_adjusted/all_available/tss5kb_correlations.tsv')}
results={label:[] for label in ['Body mass corrected','Without body mass']}
for gene,values in sorted(data.items()):
    ids=sorted(values);n=len(ids)
    if n<5:continue
    x=np.log10([float(traits[s]['MLS_years']) for s in ids])
    y=np.array([values[s] for s in ids])
    mass=np.log10([float(traits[s]['body_mass_g']) for s in ids])
    if np.ptp(x)==0 or np.ptp(y)==0:continue
    for label in results:
        X=np.column_stack([np.ones(n),x,mass] if label=='Body mass corrected' else [np.ones(n),x])
        beta,_,rank,_=np.linalg.lstsq(X,y,rcond=None)
        if rank!=X.shape[1]:continue
        df=n-rank
        se=math.sqrt(float(np.sum((y-X@beta)**2))/df*np.linalg.inv(X.T@X)[1,1])
        stat=beta[1]/se;p=float(2*t.sf(abs(stat),df));r=float(stat/math.sqrt(stat*stat+df))
        results[label].append(dict(gene=gene,model=label,n_species=n,r=r,p=p,df=df,MLS_beta=float(beta[1]),
                                   direction='Positive' if r>0 else 'Negative',species=','.join(ids)))
        if label=='Body mass corrected':
            old=prior[gene]
            assert old['species']==','.join(ids)
            for k,v in [('r',r),('p',p),('MLS_beta',float(beta[1]))]:
                assert math.isclose(v,float(old[k]),rel_tol=1e-7,abs_tol=1e-10),(gene,k)
assert {r['gene'] for r in results['Body mass corrected']}==set(prior)
summary=[]
for label,rows in results.items():
    for r,q in zip(rows,false_discovery_control([r['p'] for r in rows])):
        r['BH_FDR']=float(q)
        if label=='Body mass corrected':assert math.isclose(float(q),float(prior[r['gene']]['q']),rel_tol=1e-7,abs_tol=1e-10)
    rows.sort(key=lambda r:(r['BH_FDR'],r['p'],r['gene']))
    name='mass_corrected' if label=='Body mass corrected' else 'without_mass'
    for suffix,subset in [('all_genes',rows),('nominal_p05',[r for r in rows if r['p']<.05]),('FDR05',[r for r in rows if r['BH_FDR']<.05])]:
        write(OUT/f'{name}_{suffix}.tsv',subset,list(rows[0]))
    s=dict(model=label,tested_genes=len(rows),nominal_positive=sum(r['p']<.05 and r['r']>0 for r in rows),
           nominal_negative=sum(r['p']<.05 and r['r']<0 for r in rows),FDR_positive=sum(r['BH_FDR']<.05 and r['r']>0 for r in rows),
           FDR_negative=sum(r['BH_FDR']<.05 and r['r']<0 for r in rows),minimum_FDR=min(r['BH_FDR'] for r in rows),
           species_min=min(r['n_species'] for r in rows),species_max=max(r['n_species'] for r in rows))
    summary.append(s);print(json.dumps(s),flush=True)
write(OUT/'summary.tsv',summary,list(summary[0]))
(OUT/'validation.json').write_text(json.dumps(dict(source_checks=checks,summary=summary,
    methods='TSS ±5 kb (10 kb total). Integrated bigWig CPM·bp across the entire window, species median, log10 activity versus log10 MLS. Highest-median locus per gene and species. Assembly-span availability >=50%, positive activity, >=5 species. Human, macaque, bowhead and ambiguous/lower-bound MLS excluded. No gene-length or phylogenetic covariate. BH across all tested genes per model. Source inputs inherited from the previously validated quantification; current bigWig metadata, locus hashes, species and trait values verified. Corrected fits and BH independently reproduced by full-model OLS.',
    input_sha256=hashlib.sha256((ACT/'analysis_input_species_medians.tsv.gz').read_bytes()).hexdigest()),indent=2))
