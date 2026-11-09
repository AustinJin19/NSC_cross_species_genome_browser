#!/usr/bin/env python3
"""Validate/reuse browser CPM activity and extract matching liver expression."""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'cross_species_visual_tool/web_app').is_dir())
OUT = Path(__file__).resolve().parent.parent
WEB = ROOT/'cross_species_visual_tool/web_app'
OLD = ROOT/'output/activity_omics_comparison_2026-10-20'
sys.path[:0] = [str(WEB), str(WEB/'.deps'), str(OLD)]
import numpy as np
import quantify_activity as qa
from server import SPECS, EXCLUDED_SAMPLES
from omics import query

ALIASES = {'Ellobius lutescens':'Bramus lutescens'}
METHODS = [m for m,_ in qa.METHODS]


def write(path, rows, fields):
    path.parent.mkdir(exist_ok=True,parents=True)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt') as handle:
        w=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def fingerprint(spec):
    sid,common,latin,color,folder,prefix=spec
    loci=qa.read_loci(sid)
    files=sorted(p for p in (ROOT/'cross_species_visual_tool'/folder).glob('*.bw')
                 if p.name.lower().startswith(prefix.lower()) and
                 p.stem.removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES)
    inputs=dict(version=1,loci_sha256=hashlib.sha256(json.dumps(loci).encode()).hexdigest(),
                methods=qa.METHODS,files=[dict(path=str(p.resolve()),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in files])
    return hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest(),len(files)


def main():
    assert not query('SIRT6')['stale'], 'Expression sources changed; rebuild index first.'
    sid_for_name={s[2]:s[0] for s in SPECS}
    values=defaultdict(list)
    seen=set()
    notes=set()
    with sqlite3.connect(f'file:{WEB}/data/liver_omics.sqlite?mode=ro',uri=True) as con:
        for assay,gene,latin,sample,log,norm,norm_on,common in con.execute(
            'SELECT assay,symbol,species,sample,log,normalized,norm_on,common FROM measurements '
            'WHERE NOT(assay="rna" AND sample="BMR_Turk_Lu")'):
            sid=sid_for_name.get(ALIASES.get(latin,latin))
            if sid is None:
                continue
            gene=gene.upper()
            key=assay,gene,sid,sample
            assert key not in seen,key
            seen.add(key)
            val=log if assay=='rna' else norm
            if val is None or not math.isfinite(val):
                continue
            if assay=='protein':
                assert norm_on=='log.quant',norm_on
            values[(assay,gene,sid)].append((sample,float(val)))
            if latin=='Acomys cahirinus' and common=='Acomys russatus':
                notes.add('Acomys RNA species field says cahirinus but common field says russatus. Follow source species field; interpret the ASM comparison with this unresolved identity mismatch.')
        sources=json.loads(con.execute("SELECT value FROM metadata WHERE key='sources'").fetchone()[0])
    expr=[]
    expr_samples=defaultdict(set)
    for (assay,gene,sid),vals in sorted(values.items()):
        expr.append(dict(assay=assay,gene=gene,species=sid,median=float(np.median([v for _,v in vals])),n_samples=len(vals)))
        expr_samples[(sid,assay)].update(s for s,_ in vals)
    write(OUT/'expression_species_medians.tsv.gz',expr,['assay','gene','species','median','n_samples'])
    coverage=[]
    manifest=[]
    selected=[s for s in SPECS if any((s[0],assay) in expr_samples for assay in ['rna','protein'])]
    for spec in SPECS:
        coverage.append(dict(species=spec[0],common=spec[1],scientific_name=spec[2],
                             RNA_samples=len(expr_samples[(spec[0],'rna')]),protein_samples=len(expr_samples[(spec[0],'protein')]),
                             included=spec in selected,reason='Matching expression available' if spec in selected else 'No matching RNA or protein in the expression index'))
    write(OUT/'species_coverage.tsv',coverage,list(coverage[0]))
    qa.CACHE=OUT/'activity_cache'
    qa.CACHE.mkdir(exist_ok=True)
    for spec in selected:
        sid=spec[0]
        fp,n_samples=fingerprint(spec)
        path=None
        for cache in [OLD/'activity_cache',qa.CACHE]:
            meta=cache/f'{sid}.json'
            data=cache/f'{sid}.tsv.gz'
            if data.exists() and meta.exists() and json.loads(meta.read_text()).get('fingerprint')==fp:
                path=data
                break
        if path is None:
            print(f'{sid}: quantifying five CPM activity regions',flush=True)
            qa.run_species(spec)
            path=qa.CACHE/f'{sid}.tsv.gz'
        else:
            print(f'{sid}: verified current activity cache',flush=True)
        manifest.append(dict(species=sid,common=spec[1],scientific_name=spec[2],ATAC_samples=n_samples,
                             activity_cache=str(path),fingerprint=fp))
    (OUT/'inputs.json').write_text(json.dumps(dict(species=manifest,expression_sources=sources,
                                                notes=sorted(notes),methods=qa.METHODS),indent=2)+'\n')
    print(json.dumps(dict(species=len(selected),RNA_species=sum(bool(expr_samples[(s[0],'rna')]) for s in selected),
                         protein_species=sum(bool(expr_samples[(s[0],'protein')]) for s in selected)),indent=2),flush=True)


if __name__=='__main__':
    main()
