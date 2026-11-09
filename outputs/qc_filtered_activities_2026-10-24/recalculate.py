"""Recalculate species medians and genome-wide MLS tests after ATAC QC filtering.

Reuse verified per-sample exact integrals, never old species medians. Recount
peaks from current BED files. Preserve prior analyses and original measurements.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
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

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
WEB=ROOT/'cross_species_visual_tool/web_app'
sys.path[:0]=[str(WEB),str(WEB/'.deps')]
import activity
from server import SPECS, EXCLUDED_SAMPLES, peak_path
from sample_qc import QC_PATH, QC_POLICY, filter_records, sample_id
from calculate_gene_activity import quantify as quantify_signal

OLD=ROOT/'outputs/activity_omics_comparison_2026-10-20'
CACHE=ROOT/'output/activity_omics_comparison_2026-10-20/activity_cache'
LIFE=json.loads((WEB/'data/species_mls.json').read_text())
OLD_META={r['species']:r for r in json.loads((OLD/'input_provenance.json').read_text())['activity']}
METHODS=[('tss2kb',2000),('tss5kb',5000),('body',None),('tss200kb',200000),('tss500kb',500000)]
EXCLUDED_SPECIES={'Human','Macaque','Balaena_mysticetus'}
FIELDS=['method','body_mass_correction','gene_length_correction','gene','n_species','r','p','BH_FDR','direction','MLS_beta','df','species']


def write(path,rows,fields):
    path.parent.mkdir(exist_ok=True,parents=True)
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)


def read_records(con,sid):
    result=[]
    for table in ['activity','activity_additions']:
        result.extend(dict(r) for r in con.execute(f'SELECT * FROM {table} WHERE species=? ORDER BY locus',(sid,)))
    return result


def traits_and_records():
    traits={};records={};inventory=[]
    with sqlite3.connect(f'file:{activity.DB}?mode=ro',uri=True) as con:
        con.row_factory=sqlite3.Row
        for sid,name,latin,color,folder,prefix in SPECS:
            life=LIFE.get(latin,{})
            year,mass=life.get('years'),life.get('body_mass_g')
            rows=filter_records(read_records(con,sid))
            samples=sorted({s for r in rows for s in r['samples']})
            assert all(sample_id(s) not in EXCLUDED_SAMPLES for s in samples)
            reason=('No samples retained after QC' if not samples else
                    'Excluded by analysis settings' if sid in EXCLUDED_SPECIES else
                    'MLS missing or lower bound' if not isinstance(year,(int,float)) or year<=0 else
                    'Body mass missing' if not isinstance(mass,(int,float)) or mass<=0 else 'Included')
            row=dict(species=sid,common_name=name,scientific_name=latin,MLS_years=year,body_mass_g=mass,
                     samples_retained=len(samples),samples='; '.join(samples),status=reason)
            inventory.append(row)
            if reason=='Included':traits[sid]=row;records[sid]=rows
    return traits,records,inventory


def fit(gene,values,traits,method,mass_correction,length_correction):
    ids=sorted(s for s,r in values.items() if r['region_eligible'] and r['median']>0)
    n=len(ids)
    if n<5:return None
    x=np.log10([traits[s]['MLS_years'] for s in ids]);y=np.log10([values[s]['median'] for s in ids])
    if np.ptp(x)==0 or np.ptp(y)==0:return None
    columns=[np.ones(n)]
    if mass_correction:columns.append(np.log10([traits[s]['body_mass_g'] for s in ids]))
    if length_correction:columns.append(np.log10([values[s]['length_kb'] for s in ids]))
    C=np.column_stack(columns)
    if np.linalg.matrix_rank(C)<C.shape[1]:return None
    Q,_=np.linalg.qr(C,mode='reduced');rx=x-Q@(Q.T@x);ry=y-Q@(Q.T@y)
    xx,yy=rx@rx,ry@ry
    if xx<=1e-12*np.sum((x-x.mean())**2) or yy<=1e-12*np.sum((y-y.mean())**2):return None
    r=float(np.clip(rx@ry/math.sqrt(xx*yy),-1,1));df=n-C.shape[1]-1
    if df<1:return None
    p=float(2*t.sf(abs(r)*math.sqrt(df/max(1-r*r,1e-300)),df))
    X=np.column_stack([C,x]);coef=np.linalg.lstsq(X,y,rcond=None)[0]
    error=y-X@coef;covariance=(error@error)/df*np.linalg.inv(X.T@X)
    ols_p=float(2*t.sf(abs(coef[-1]/math.sqrt(covariance[-1,-1])),df))
    assert abs(p-ols_p)<1e-8,(gene,method,p,ols_p)
    return dict(method=method,body_mass_correction=mass_correction,gene_length_correction=length_correction,
                gene=gene,n_species=n,r=r,p=p,direction='Positive' if r>0 else 'Negative',
                MLS_beta=float(coef[-1]),df=df,species=','.join(ids),ols_p_difference=abs(p-ols_p))


def export_fits(family,genes,traits):
    summary=[];all_sig=[];max_difference=0.
    for method,gene_values in genes.items():
        for length in ([True,False] if method=='body' else [False]):
            for mass in [True,False]:
                rows=[r for gene,values in sorted(gene_values.items())
                      if (r:=fit(gene,values,traits,method,mass,length)) is not None]
                for row,q in zip(rows,false_discovery_control([r['p'] for r in rows])):row['BH_FDR']=float(q)
                rows.sort(key=lambda r:(r['p'],r['gene']))
                name=f'{method}_mass{int(mass)}_length{int(length)}'
                dest=OUT/family/name
                write(dest/'all_tested.tsv.gz',rows,FIELDS)
                for cutoff,key in [('FDR05','BH_FDR'),('nominal_p05','p')]:
                    selected=[r for r in rows if r[key]<.05]
                    write(dest/f'{cutoff}.tsv',selected,FIELDS)
                    for sign in ['Positive','Negative']:
                        sub=[r for r in selected if r['direction']==sign]
                        write(dest/f'{cutoff}_{sign.lower()}.tsv',sub,FIELDS)
                sig=[r for r in rows if r['BH_FDR']<.05];all_sig.extend(sig)
                minimum=min(r['BH_FDR'] for r in rows)
                summary.append(dict(family=family,method=method,body_mass_correction=mass,gene_length_correction=length,
                    tested=len(rows),FDR_positive=sum(r['r']>0 for r in sig),FDR_negative=sum(r['r']<0 for r in sig),
                    nominal_positive=sum(r['p']<.05 and r['r']>0 for r in rows),
                    nominal_negative=sum(r['p']<.05 and r['r']<0 for r in rows),minimum_FDR=minimum))
                max_difference=max(max_difference,max(r['ols_p_difference'] for r in rows))
                print(json.dumps(summary[-1]),flush=True)
    write(OUT/f'{family}_summary.tsv',summary,list(summary[0]))
    write(OUT/f'{family}_significant_genes.tsv',all_sig,FIELDS)
    return summary,max_difference


def run_activity(traits,records):
    genes={m:defaultdict(dict) for m,k in METHODS};checks=[]
    fields=['method','gene','species','locus','median','length_kb','available_fraction','region_eligible','n_samples','samples']
    for sid,trait in traits.items():
        meta=OLD_META[sid]
        # Cached sample arrays are aligned to all old loci, even for excluded samples.
        with sqlite3.connect(f'file:{activity.DB}?mode=ro',uri=True) as con:
            con.row_factory=sqlite3.Row
            loci=sorted({(r['symbol_key'],r['locus']) for r in read_records(con,sid)},key=lambda r:r[1])
        assert hashlib.sha256(json.dumps(loci).encode()).hexdigest()==meta['inputs']['loci_sha256'],sid
        assert meta['inputs']['methods']==[list(m) for m in METHODS]
        matrices=[];fractions=[];sample_names=[]
        n=len(loci);requests=[];lengths=[];regions=[]
        for method,flank in METHODS:
            for gene,locus in loci:
                chrom,span,strand=locus.rsplit('|',1)[1].rsplit(':',2);a,b=map(int,span.split('-'))
                tss=a if strand=='+' else b-1
                start,end=((a-2000,b) if strand=='+' else (a,b+2000)) if flank is None else (tss-flank,tss+flank)
                regions.append(dict(chrom=chrom,region_start=start,region_end=end))
                requests.append(end-start);lengths.append((b-a)/1000)
        requests=np.asarray(requests)
        allowed={sample_id(s) for r in records[sid] for s in r['samples']}
        source=[]
        for info in meta['inputs']['files']:
            path=Path(info['path']);st=path.stat()
            assert (st.st_size,st.st_mtime_ns)==(info['size'],info['mtime_ns']),str(path)
            if sample_id(path) not in allowed:continue
            assert sample_id(path) not in EXCLUDED_SAMPLES
            with np.load(CACHE/f'{sid}__{path.stem}.npz') as cache:
                assert str(cache['fingerprint'])==meta['fingerprint']
                signal,width=cache['signal'].copy(),cache['widths'].copy()
            assert len(signal)==len(regions)
            # Check retained cached sample values against current exact bigWig integration.
            ix=np.linspace(0,len(regions)-1,15,dtype=int)
            direct,_,dw,_=quantify_signal(path,[regions[i] for i in ix])
            assert np.allclose(signal[ix],direct,rtol=2e-7,atol=1e-5,equal_nan=True),(sid,path.stem)
            assert np.array_equal(width[ix],dw)
            matrices.append(signal);fractions.append(np.where(np.isfinite(signal),width/requests,np.nan))
            sample_names.append(path.stem);source.append(info)
        assert set(map(sample_id,sample_names))==allowed
        values=np.asarray(matrices)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',category=RuntimeWarning)
            medians=np.nanmedian(values,axis=0);available=np.nanmin(np.asarray(fractions),axis=0)
        n_samples=np.isfinite(values).sum(axis=0)
        selected_rows=[]
        for k,(method,_) in enumerate(METHODS):
            selected={}
            for i,(gene,locus) in enumerate(loci):
                j=k*n+i;v=medians[j]
                if not math.isfinite(v):continue
                if gene not in selected or v>selected[gene]['median']:
                    selected[gene]=dict(method=method,gene=gene,species=sid,locus=locus,median=float(v),
                        length_kb=lengths[j],available_fraction=float(available[j]),region_eligible=bool(available[j]>=.5),
                        n_samples=int(n_samples[j]),samples=json.dumps({s:float(values[z,j]) if np.isfinite(values[z,j]) else None for z,s in enumerate(sample_names)}))
            for gene,row in selected.items():genes[method][gene][sid]=row
            selected_rows.extend(selected.values())
        write(OUT/'activity_selected_loci'/f'{sid}.tsv.gz',selected_rows,fields)
        checks.append(dict(species=sid,samples=sample_names,source_files=source,loci=len(loci),sample_integrals_checked=15*len(sample_names)))
        print(f'Activity: {sid}: {len(sample_names)} samples; recomputed all five species-median methods',flush=True)
    summary,err=export_fits('gene_activity',genes,traits)
    return dict(summary=summary,source_checks=checks,max_OLS_p_difference=err)


def run_regulatory(traits,records):
    path=ROOT/'outputs/regulatory_windows_2026-10-22/analyze_windows.py'
    spec=importlib.util.spec_from_file_location('regulatory_reference',path)
    ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
    genes={f'peaks{kb}kb':defaultdict(dict) for kb in ref.WINDOWS}
    folders={s[0]:ROOT/'cross_species_visual_tool'/s[4] for s in SPECS}
    checks=[]
    for sid,trait in traits.items():
        rows=sorted(records[sid],key=lambda r:r['locus'])
        medians,available,n_samples,n_check=ref.quantify(rows,sid,folders[sid])
        sources=[]
        for sample in sorted({s for r in rows for s in r['samples']}):
            bw=folders[sid]/(sample+'.bw')
            assert sample_id(sample) not in EXCLUDED_SAMPLES
            for path in [bw,peak_path(bw)]:
                st=path.stat() if path.exists() else None
                sources.append(dict(path=str(path),size=st.st_size if st else None,mtime_ns=st.st_mtime_ns if st else None))
        for k,kb in enumerate(ref.WINDOWS):
            selected={};method=f'peaks{kb}kb'
            for i,row in enumerate(rows):
                v=medians[k,i];gene=row['symbol_key']
                if not math.isfinite(v):continue
                if gene not in selected or v>selected[gene]['median']:
                    selected[gene]=dict(method=method,gene=gene,species=sid,locus=row['locus'],median=float(v),
                        available_fraction=float(available[k,i]),region_eligible=bool(available[k,i]>=.5),n_samples=int(n_samples[k,i]))
            for gene,row in selected.items():genes[method][gene][sid]=row
            write(OUT/'regulatory_selected_loci'/f'{sid}_{kb}kb.tsv.gz',selected.values(),
                  ['method','gene','species','locus','median','available_fraction','region_eligible','n_samples'])
        checks.append(dict(species=sid,source_files=sources,browser_scorer_checks=n_check))
        print(f'Regulatory: {sid}: recounted all six peak windows',flush=True)
    summary,err=export_fits('regulatory_activity',genes,traits)
    return dict(summary=summary,source_checks=checks,max_OLS_p_difference=err)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['activity','regulatory']);args=parser.parse_args()
    started=time.time();traits,records,inventory=traits_and_records()
    write(OUT/f'{args.mode}_species.tsv',inventory,list(inventory[0]))
    report=run_activity(traits,records) if args.mode=='activity' else run_regulatory(traits,records)
    report.update(qc_policy=QC_POLICY,qc_sha256=hashlib.sha256(QC_PATH.read_bytes()).hexdigest(),
        traits_sha256=hashlib.sha256((WEB/'data/species_mls.json').read_bytes()).hexdigest(),
        activity_database_mtime_ns=activity.DB.stat().st_mtime_ns,elapsed_seconds=time.time()-started,
        method='log10(species median) versus log10(MLS), mass on/off; body also length on/off; >=5 species; >=50% available span; highest median locus before eligibility; BH separately within each method/model',
        excluded_species=sorted(EXCLUDED_SPECIES),eligible_species=len(traits))
    (OUT/f'{args.mode}_manifest.json').write_text(json.dumps(report,indent=2))
    print('Completed',args.mode,'in',round(time.time()-started,1),'seconds',flush=True)
