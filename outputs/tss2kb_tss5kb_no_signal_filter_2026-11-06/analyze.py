"""Genome-wide TSS activity/MLS scan without a positive-signal coverage filter.

Reproducible run: python3 outputs/tss2kb_tss5kb_no_signal_filter_2026-11-06/analyze.py
Uses exact bigWig interval integrals, current catalog samples and name mapping.
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
APP = ROOT/'cross_species_visual_tool/web_app'
sys.path[:0] = [str(APP), str(APP/'.deps')]
import numpy as np
import pyBigWig
from scipy.stats import t
import activity
import gene_names
import server
from sample_qc import EXCLUDED_SAMPLES, QC_POLICY, sample_id
from signal_coverage import positive_bases

METHODS = [('tss2kb', 2000), ('tss5kb', 5000)]
EXCLUDED = {'Human', 'Balaena_mysticetus'}
MIN_SPECIES = 5
LIFE = json.loads((APP/'data/species_mls.json').read_text())
CACHE = OUT.parent/'tss2kb_tss5kb_positive_signal_2026-11-06/quantification'
FIELDS = ['gene','window_kb','body_mass_correction','n_species','direction','r','p','BH_FDR',
          'MLS_beta','MLS_beta_SE','mass_beta','df','species']


def write_tsv(path, rows, fields):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def fingerprint(path):
    s = path.stat()
    return dict(path=str(path.resolve()), size=s.st_size, mtime_ns=s.st_mtime_ns)


def read_loci(sid):
    with sqlite3.connect(f'file:{activity.DB}?mode=ro', uri=True) as con:
        loci = {r for table in ('activity', 'activity_additions') for r in con.execute(
            f'SELECT symbol_key,locus FROM {table} WHERE species=?', (sid,))}
    return sorted(loci, key=lambda r:r[1])


def regions(loci):
    result = []
    for method, flank in METHODS:
        for symbol, locus in loci:
            chrom, span, strand = locus.rsplit('|', 1)[1].rsplit(':', 2)
            a, b = map(int, span.split('-'))
            tss = a if strand == '+' else b-1
            result.append((chrom, tss-flank, tss+flank))
    return result


def integrate_block(records, lo, hi, qs, qe):
    """Piecewise integrals for total CPM and finite positive base count."""
    if not records:
        return np.zeros(len(qs)), np.zeros(len(qs))
    arr = np.asarray(records, dtype=float)
    arr = arr[np.isfinite(arr[:, 2])]
    if not len(arr): return np.zeros(len(qs)), np.zeros(len(qs))
    a, b, v = arr.T
    a = np.maximum(a, lo); b = np.minimum(b, hi)
    assert np.all(a[1:] >= b[:-1])
    positive = (v > 0).astype(float)
    signal = np.r_[0., np.cumsum((b-a)*v)]
    coverage = np.r_[0., np.cumsum((b-a)*positive)]
    def at(x, cumulative, values):
        j = np.searchsorted(b, x, side='right')
        total = cumulative[j].copy(); mask = j < len(a)
        k = j[mask]
        total[mask] += np.maximum(0, x[mask]-a[k])*values[k]
        return total
    return at(qe, signal, v)-at(qs, signal, v), at(qe, coverage, positive)-at(qs, coverage, positive)


def quantify(path, requests):
    count = len(requests)
    scores = np.full(count, np.nan); positives = np.full(count, np.nan)
    widths = np.zeros(count, dtype=int)
    with pyBigWig.open(str(path)) as bw:
        sizes = bw.chroms(); grouped = defaultdict(list)
        for i, (c, a, b) in enumerate(requests):
            chrom = server.resolve_chrom(c, sizes)
            if chrom is not None: grouped[chrom].append(i)
        for chrom, indexes in grouped.items():
            indexes = np.asarray(indexes)
            starts = np.array([max(0, requests[i][1]) for i in indexes])
            ends = np.array([min(sizes[chrom], requests[i][2]) for i in indexes])
            valid = ends > starts
            ids, qs, qe = indexes[valid], starts[valid], ends[valid]
            if not len(ids): continue
            signal = np.zeros(len(ids)); positive = np.zeros(len(ids))
            block = 1_000_000
            occupied = sorted({j for a, b in zip(qs, qe) for j in range(int(a)//block, (int(b)-1)//block+1)})
            for j in occupied:
                lo = j*block; hi = min(lo+block, sizes[chrom])
                hits = np.flatnonzero((qs < hi) & (qe > lo))
                s, p = integrate_block(bw.intervals(chrom, lo, hi), lo, hi,
                                       np.maximum(qs[hits], lo), np.minimum(qe[hits], hi))
                signal[hits] += s; positive[hits] += p
            # Floating subtraction can produce tiny round-off around zero.
            signal[np.abs(signal) < 1e-9] = 0.
            scores[ids] = signal; positives[ids] = positive; widths[ids] = qe-qs
    return scores, positives, widths


def run_species(spec):
    sid, common, latin, color, folder, prefix = spec; begin = time.time()
    loci = read_loci(sid); requests = regions(loci)
    files = sorted(p for p in (server.ROOT/folder).glob('*.bw')
                   if p.name.lower().startswith(prefix.lower()) and sample_id(p) not in EXCLUDED_SAMPLES)
    inputs = dict(version=1, loci_sha256=hashlib.sha256(json.dumps(loci).encode()).hexdigest(),
                  windows=METHODS, files=[fingerprint(p) for p in files])
    fp = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    matrices = []; coverage = []; fractions = []; checks = []
    requested = np.array([b-a for c,a,b in requests])
    for path in files:
        dest = CACHE/f'{sid}__{path.stem}.npz'
        if dest.exists():
            with np.load(dest) as z:
                if str(z['fingerprint']) != fp: raise ValueError(f'Stale cache: {dest}')
                scores, positive, width = z['scores'], z['positive_bp'], z['widths']
        else:
            scores, positive, width = quantify(path, requests)
            np.savez_compressed(dest, fingerprint=fp, scores=scores, positive_bp=positive, widths=width)
        # Independent direct exact bigWig queries, including TERT and SIRT6.
        rng = np.random.default_rng(1106)
        selected = set(rng.choice(len(requests), min(20, len(requests)), replace=False))
        selected.update(i + k*len(loci) for i,(s,l) in enumerate(loci) if s in ('TERT','SIRT6') for k in range(2))
        with pyBigWig.open(str(path)) as bw:
            sizes = bw.chroms(); max_error = 0.
            for i in sorted(selected):
                c, a, b = requests[i]; c = server.resolve_chrom(c, sizes)
                if c is None or min(b,sizes[c]) <= max(0,a):
                    assert np.isnan(scores[i]); continue
                direct = bw.stats(c,max(0,a),min(b,sizes[c]),type='sum',exact=True)[0] or 0.
                exact_positive = positive_bases(str(path),path.stat().st_size,path.stat().st_mtime_ns,c,a,b)
                assert np.isclose(scores[i],direct,rtol=2e-7,atol=1e-5), (sid,path.stem,i,scores[i],direct)
                assert positive[i] == exact_positive, (sid,path.stem,i,positive[i],exact_positive)
                max_error = max(max_error,abs(scores[i]-direct))
        checks.append(dict(sample=path.stem,regions=len(selected),maximum_signal_difference=max_error))
        matrices.append(scores); coverage.append(positive/requested)
        fractions.append(np.where(np.isfinite(scores),width/requested,np.nan))
        print(f'{sid}: {path.stem}, {len(loci):,} loci x 2 windows, {time.time()-begin:.1f}s',flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        medians = np.nanmedian(matrices,axis=0)
        positive_fraction = np.nanmedian(coverage,axis=0)
        available_fraction = np.nanmin(fractions,axis=0)
    np.savez_compressed(CACHE/f'{sid}.npz',median=medians,positive_fraction=positive_fraction,
                        available_fraction=available_fraction,n_samples=np.isfinite(matrices).sum(axis=0))
    (CACHE/f'{sid}.json').write_text(json.dumps(dict(inputs=inputs,fingerprint=fp,loci=loci,
                            samples=[p.stem for p in files],checks=checks),indent=2))
    print(f'{sid}: COMPLETE {time.time()-begin:.1f}s',flush=True)
    return sid


def fit(gene, observations, traits, method, mass):
    ids = sorted(s for s,r in observations.items() if r['eligible'] and r['median_activity']>0
                 and (not mass or traits[s]['body_mass_g']>0))
    n = len(ids)
    if n < MIN_SPECIES: return None
    x = np.log10([traits[s]['MLS_years'] for s in ids])
    y = np.log10([observations[s]['median_activity'] for s in ids])
    if np.ptp(x)==0 or np.ptp(y)==0:return None
    C = np.column_stack([np.ones(n)] + ([np.log10([traits[s]['body_mass_g'] for s in ids])] if mass else []))
    if np.linalg.matrix_rank(C)<C.shape[1]:return None
    Q,_ = np.linalg.qr(C,mode='reduced'); rx=x-Q@(Q.T@x); ry=y-Q@(Q.T@y)
    xx, yy = rx@rx, ry@ry
    if xx<=1e-12*np.sum((x-x.mean())**2) or yy<=1e-12*np.sum((y-y.mean())**2):return None
    r = float(np.clip((rx@ry)/math.sqrt(xx*yy),-1,1)); df=n-C.shape[1]-1
    p = float(2*t.sf(abs(r)*math.sqrt(df/max(1-r*r,1e-300)),df))
    X = np.column_stack([C,x]); coef=np.linalg.lstsq(X,y,rcond=None)[0]
    error=y-X@coef; covariance=(error@error)/df*np.linalg.inv(X.T@X)
    se=math.sqrt(covariance[-1,-1]); coef_p=float(2*t.sf(abs(coef[-1]/se),df))
    assert np.isclose(p,coef_p,atol=1e-8), (gene,p,coef_p)
    return dict(gene=gene,window_kb=int(method[3:-2]),body_mass_correction='Yes' if mass else 'No',
                n_species=n,direction='Positive' if r>0 else 'Negative',r=r,p=p,
                MLS_beta=float(coef[-1]),MLS_beta_SE=se,mass_beta=float(coef[1]) if mass else None,
                df=df,species='; '.join(ids))


def bh(values):
    p = np.asarray(values); order=np.argsort(p); ranked=p[order]*len(p)/np.arange(1,len(p)+1)
    q=np.minimum.accumulate(ranked[::-1])[::-1].clip(0,1); out=np.empty(len(p)); out[order]=q
    return out


def analyze(specs, inventory):
    resolver=gene_names.resolver(); sources={}; symbol_counts=defaultdict(set)
    for spec in specs:
        sid=spec[0]; meta=json.loads((CACHE/f'{sid}.json').read_text()); source=dict(np.load(CACHE/f'{sid}.npz'))
        source.update(meta); by_symbol=defaultdict(list)
        for i,(s,l) in enumerate(meta['loci']): by_symbol[s].append(i); symbol_counts[resolver.canonical(s)].add(sid)
        source['by_symbol']=by_symbol; sources[sid]=source
    candidates=sorted(symbol_counts)
    # Resolve against the same original symbols as the browser, without merging ambiguous paralogs.
    mapping={}; evidence=[]
    for gene in candidates:
        matches={}
        for sid,source in sources.items():
            m=resolver.resolve(gene,sid,source['by_symbol'])
            if m['symbol'] is not None:
                matches[sid]=m
        if len(matches)>=MIN_SPECIES:
            mapping[gene]=matches
    print(f'{len(mapping):,} canonical candidate genes with annotation matches in >= {MIN_SPECIES} species',flush=True)
    traits={r['species']:r for r in inventory if r['status']=='Included'}
    genes={method:defaultdict(dict) for method,_ in METHODS}
    details=[]
    for gene,matches in mapping.items():
        for sid,match in matches.items():
            source=sources[sid]; n=len(source['loci'])
            for k,(method,_) in enumerate(METHODS):
                indexes=[i+k*n for i in source['by_symbol'][match['symbol']]]
                finite=[i for i in indexes if np.isfinite(source['median'][i])]
                if not finite:continue
                # Select highest median BEFORE the genomic-span filter, lexical tie break.
                i=min(finite,key=lambda j:(-source['median'][j],source['loci'][j%n][1]))
                value=float(source['median'][i]); pos=float(source['positive_fraction'][i]); avail=float(source['available_fraction'][i])
                eligible=bool(np.isfinite(avail) and avail>=.5)
                row=dict(gene=gene,method=method,species=sid,locus=source['loci'][i%n][1],
                         native_symbol=match['symbol'],match_method=match['method'],median_activity=value,
                         positive_signal_fraction=pos,available_fraction=avail,eligible=eligible,
                         n_samples=int(source['n_samples'][i]))
                genes[method][gene][sid]=row; details.append(row)
    write_tsv(OUT/'selected_loci_and_coverage.tsv.gz',details,list(details[0]))
    summary=[]; all_results=[]
    for method,_ in METHODS:
        for mass in (False,True):
            rows=[r for g,v in genes[method].items() if (r:=fit(g,v,traits,method,mass)) is not None]
            for r,q in zip(rows,bh([r['p'] for r in rows])):r['BH_FDR']=float(q)
            rows.sort(key=lambda r:(r['p'],r['gene'])); all_results.extend(rows)
            stem=method+('_with_body_mass' if mass else '_without_body_mass')
            write_tsv(OUT/(stem+'_all_tested.tsv.gz'),rows,FIELDS)
            for label,key in [('FDR05','BH_FDR'),('nominal_p05','p')]:
                selected=[r for r in rows if r[key]<.05]
                write_tsv(OUT/(stem+'_'+label+'.tsv'),selected,FIELDS)
                for sign in ('Positive','Negative'):
                    names=[r['gene'] for r in selected if r['direction']==sign]
                    (OUT/(stem+'_'+label+'_'+sign.lower()+'.txt')).write_text('\n'.join(names)+ ('\n' if names else ''))
            summary.append(dict(window_kb=int(method[3:-2]),body_mass_correction='Yes' if mass else 'No',
                                genes_tested=len(rows),FDR_positive=sum(r['r']>0 and r['BH_FDR']<.05 for r in rows),
                                FDR_negative=sum(r['r']<0 and r['BH_FDR']<.05 for r in rows),
                                nominal_positive=sum(r['r']>0 and r['p']<.05 for r in rows),
                                nominal_negative=sum(r['r']<0 and r['p']<.05 for r in rows),
                                min_species=min(r['n_species'] for r in rows),max_species=max(r['n_species'] for r in rows)))
            print(json.dumps(summary[-1]),flush=True)
    write_tsv(OUT/'summary.tsv',summary,list(summary[0]))
    (OUT/'results.json').write_text(json.dumps(dict(summary=summary,results=all_results,inventory=inventory)))
    for mass in ('No','Yes'):
        selected=[r for r in all_results if r['body_mass_correction']==mass and r['BH_FDR']<.05]
        write_tsv(OUT/('Significant_genes_'+('with' if mass=='Yes' else 'without')+'_body_mass.tsv'),selected,FIELDS)
    report=dict(date='2026-11-06',model='log10(species median integrated CPM activity) ~ log10(MLS) + optional log10(body mass in g)',
                regression='Ordinary least squares; two-sided Pearson/partial-correlation t test; no phylogenetic correction',
                minimum_species=MIN_SPECIES,FDR='Benjamini-Hochberg separately within each window and model over all tested genes',
                positive_signal_filter='Disabled: positive-base coverage is not an eligibility criterion; nonpositive median activity remains excluded for log10',
                available_region_filter='Minimum available scaffold fraction across measured samples >= 0.5',
                locus_selection='Highest species median before coverage filtering; lexical locus tie break',
                naming='Current offline browser HGNC alias / HCOP / projected-symbol resolver; ambiguous mappings not pooled',
                excluded_species=sorted(EXCLUDED),qc_policy=QC_POLICY,summary=summary,
                quantification_cache=str(CACHE),
                cache_validation='Current annotation loci, windows and every sample bigWig path/size/mtime fingerprint verified; species aggregates recomputed from sample caches and checked',
                inputs=[fingerprint(activity.DB),fingerprint(gene_names.DATA),fingerprint(APP/'data/species_mls.json')])
    (OUT/'analysis_manifest.json').write_text(json.dumps(report,indent=2))
    print('Analysis complete.',flush=True)


def verify_cache(spec):
    """Reuse immutable measurements only after validating current input sources."""
    sid, common, latin, color, folder, prefix = spec
    loci = read_loci(sid)
    files = sorted(p for p in (server.ROOT/folder).glob('*.bw')
                   if p.name.lower().startswith(prefix.lower()) and sample_id(p) not in EXCLUDED_SAMPLES)
    inputs = dict(version=1, loci_sha256=hashlib.sha256(json.dumps(loci).encode()).hexdigest(),
                  windows=METHODS, files=[fingerprint(p) for p in files])
    expected = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
    meta = json.loads((CACHE/f'{sid}.json').read_text())
    assert meta['fingerprint'] == expected, f'Stale source cache: {sid}'
    assert meta['loci'] == [list(row) for row in loci]
    sample_scores, available = [], []
    requested = np.array([b-a for c,a,b in regions(loci)])
    for path in files:
        with np.load(CACHE/f'{sid}__{path.stem}.npz') as z:
            assert str(z['fingerprint']) == expected, f'Stale sample cache: {path}'
            sample_scores.append(z['scores'])
            available.append(np.where(np.isfinite(z['scores']),z['widths']/requested,np.nan))
    with np.load(CACHE/f'{sid}.npz') as z, warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        np.testing.assert_allclose(z['median'],np.nanmedian(sample_scores,axis=0),equal_nan=True)
        np.testing.assert_allclose(z['available_fraction'],np.nanmin(available,axis=0),equal_nan=True)
        np.testing.assert_array_equal(z['n_samples'],np.isfinite(sample_scores).sum(axis=0))
    print(f'{sid}: verified {len(files)} unchanged sample sources and {len(loci):,} loci',flush=True)


def main():
    inventory=[]; specs=[]
    for spec in server.SPECS:
        sid,common,latin,*_=spec; life=LIFE.get(latin,{})
        years=life.get('years'); mass=life.get('body_mass_g')
        status='Excluded by settings' if sid in EXCLUDED else 'MLS missing or lower bound' if not years or years<=0 else 'Included'
        inventory.append(dict(species=sid,common_name=common,scientific_name=latin,MLS_years=years,body_mass_g=mass,status=status))
        if status=='Included':specs.append(spec)
    write_tsv(OUT/'species_inventory.tsv',inventory,list(inventory[0]))
    for spec in specs:
        verify_cache(spec)
    analyze(specs,inventory)


if __name__=='__main__':main()
