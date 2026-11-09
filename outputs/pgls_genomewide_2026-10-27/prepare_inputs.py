"""Freeze browser-compatible species medians for the PGLS genome-wide scan."""
import csv,gzip,hashlib,json,math,sqlite3,sys,time,warnings
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[1]; APP=ROOT/'cross_species_visual_tool/web_app'
sys.path[:0]=[str(APP),str(APP/'.deps')]
from activity import DB
from server import SPECS,peak_path
from sample_qc import filter_records,QC_POLICY,sample_id
from calculate_gene_activity import quantify

ACT=ROOT/'output/activity_omics_comparison_2026-10-20/activity_cache'
REG=ROOT/'outputs/regulatory_windows_2026-10-22'
LIFE=json.loads((APP/'data/species_mls.json').read_text())
EXCLUDED={'Human','Macaque','Balaena_mysticetus'}
METHODS=['tss2kb','tss5kb','body','tss100kb','tss200kb','tss500kb']+[f'peaks{k}kb' for k in (20,50,100,200,500,1000)]

def records_for(sid):
    with sqlite3.connect(f'file:{DB}?mode=ro',uri=True) as con:
        con.row_factory=sqlite3.Row
        raw=[dict(r) for table in ['activity','activity_additions'] for r in con.execute(f'SELECT * FROM {table} WHERE species=? ORDER BY locus',(sid,))]
    return filter_records(raw)

def file_check(info):
    p=Path(info['path']);s=p.stat() if p.exists() else None
    assert (s.st_size if s else None,s.st_mtime_ns if s else None)==(info['size'],info['mtime_ns']),f'Stale source: {p}'

def quantify_100kb(spec):
    sid,common,latin,color,folder,prefix=spec; started=time.time()
    path=OUT/'tss100kb_cache'/f'{sid}.tsv.gz'; meta_path=path.with_suffix('.json')
    rows=records_for(sid);rows.sort(key=lambda r:r['locus'])
    samples=sorted({s for r in rows for s in r['samples']})
    files=[ROOT/'cross_species_visual_tool'/folder/(s+'.bw') for s in samples]
    metadata=dict(loci_sha256=hashlib.sha256(json.dumps([(r['symbol_key'],r['locus']) for r in rows]).encode()).hexdigest(),
                  files=[dict(path=str(p.resolve()),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in files])
    if path.exists() and meta_path.exists() and json.loads(meta_path.read_text())['inputs']==metadata:return sid
    requests=[]
    for r in rows:
        c,span,strand=r['locus'].rsplit('|',1)[1].rsplit(':',2);a,b=map(int,span.split('-'))
        tss=a if strand=='+' else b-1
        requests.append(dict(chrom=c,region_start=tss-100000,region_end=tss+100000))
    arrays=[];widths=[];checks=[]
    for sample,p in zip(samples,files):
        v,_,w,_=quantify(p,requests,chunked=True)
        present=np.array([sample in r['samples'] for r in rows]);v[~present]=np.nan
        ix=np.linspace(0,len(rows)-1,15,dtype=int)
        direct,_,dw,_=quantify(p,[requests[i] for i in ix])
        assert np.allclose(v[ix],direct,rtol=2e-7,atol=1e-5,equal_nan=True)
        assert np.array_equal(w[ix],dw)
        arrays.append(v);widths.append(np.where(np.isfinite(v),w/200000,np.nan))
        checks.append(dict(sample=sample,intervals=len(ix)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        median=np.nanmedian(arrays,axis=0);available=np.nanmin(widths,axis=0)
    chosen={}
    for i,r in enumerate(rows):
        gene=r['symbol_key'];v=median[i]
        if np.isfinite(v) and (gene not in chosen or v>median[chosen[gene]]):chosen[gene]=i
    with gzip.open(path,'wt') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(['method','gene','species','locus','median_activity','available_fraction','region_eligible'])
        for gene,i in sorted(chosen.items()):writer.writerow(['tss100kb',gene,sid,rows[i]['locus'],median[i],available[i],int(available[i]>=.5)])
    meta_path.write_text(json.dumps(dict(inputs=metadata,checks=checks),indent=2))
    print(sid,'100 kb quantified',round(time.time()-started,1),'seconds',flush=True)
    return sid

def main():
    started=time.time();(OUT/'tss100kb_cache').mkdir(exist_ok=True);(OUT/'inputs').mkdir(exist_ok=True)
    tree=json.loads((APP/'data/phylogenetic_covariance.json').read_text())
    selected=[];species=[];inventory=[];counts=Counter();validation=[]
    reg_meta=json.loads((REG/'manifest.json').read_text())
    assert reg_meta['database_mtime_ns']==DB.stat().st_mtime_ns
    assert reg_meta['traits_sha256']==hashlib.sha256((APP/'data/species_mls.json').read_bytes()).hexdigest()
    for info in reg_meta['source_files']:file_check(info)
    for spec in SPECS:
        sid,common,latin,color,folder,prefix=spec;life=LIFE.get(latin,{})
        years,mass=life.get('years'),life.get('body_mass_g')
        reason=('Excluded by settings' if sid in EXCLUDED else 'MLS missing or lower bound' if not years or years<=0 else 'Body mass unavailable' if not mass or mass<=0 else 'Absent from TimeTree' if latin not in tree['species'] else 'Included')
        inventory.append(dict(species=sid,common_name=common,scientific_name=latin,MLS_years=years,body_mass_g=mass,status=reason))
        if reason!='Included':continue
        rows=records_for(sid);allowed={s for r in rows for s in r['samples']}
        meta=json.loads((ACT/f'{sid}.json').read_text())
        loci=sorted({(r['symbol_key'],r['locus']) for r in rows},key=lambda r:r[1])
        assert hashlib.sha256(json.dumps(loci).encode()).hexdigest()==meta['inputs']['loci_sha256'],sid
        assert allowed=={Path(v['path']).stem for v in meta['inputs']['files']},(sid,allowed)
        for info in meta['inputs']['files']:file_check(info)
        # Regulatory caches have no active-TSS-filter omissions and unchanged inputs.
        reg_samples={v['sample'] for v in reg_meta['source_files'] if v['species']==sid and v['size'] is not None}
        assert reg_samples==allowed,(sid,reg_samples,allowed)
        selected.append(spec);species.append(inventory[-1]);counts.update({r['symbol_key'] for r in rows})
        validation.append(dict(species=sid,files_verified=len(meta['inputs']['files']),samples=sorted(allowed),loci=len(loci),activity_cache_fingerprint=meta['fingerprint']))
    (OUT/'species_inventory.json').write_text(json.dumps(inventory,indent=2))
    (OUT/'species.json').write_text(json.dumps(species,indent=2))
    # A gene present in fewer than five species cannot pass any scan's n filter.
    genes=sorted(g for g,n in counts.items() if n>=5);index={g:i for i,g in enumerate(genes)}
    ix=[tree['species'].index(s['scientific_name']) for s in species]
    C=np.asarray(tree['covariance'])[np.ix_(ix,ix)]
    base=dict(genes=np.asarray(genes),species=np.asarray([s['species'] for s in species]),
              scientific_names=np.asarray([s['scientific_name'] for s in species]),
              x=np.log10([s['MLS_years'] for s in species]),mass=np.log10([s['body_mass_g'] for s in species]),covariance=C)
    shape=(len(genes),len(species)); values={m:np.full(shape,np.nan) for m in METHODS};length=np.full(shape,np.nan)
    eligible_counts={m:0 for m in METHODS}
    def consume(path,j,regulatory=False):
        with gzip.open(path,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                gene=r['gene']
                if gene not in index:continue
                method=r.get('method') or 'peaks'+path.name.split('_')[-1].replace('.tsv.gz','')
                if method not in values:continue
                raw=float(r['median'] if regulatory else r['median_activity']);available=float(r['available_fraction'])
                if available<.5 or not math.isfinite(available) or raw<=0 or not math.isfinite(raw):continue
                i=index[gene];values[method][i,j]=math.log10(raw);eligible_counts[method]+=1
                if method=='body':length[i,j]=float(r['length_kb'])
    for j,s in enumerate(species):
        sid=s['species'];consume(ACT/f'{sid}.tsv.gz',j)
        for kb in (20,50,100,200,500,1000):consume(REG/'selected_loci'/f'{sid}_{kb}kb.tsv.gz',j,True)
    # Save available inputs immediately so their fits can run during quantification.
    for method in METHODS:
        if method=='tss100kb':continue
        np.savez_compressed(OUT/'inputs'/f'{method}.npz',**base,y=values[method],length_kb=length)
    (OUT/'input_validation.json').write_text(json.dumps(dict(species=validation,regulatory_source_files_verified=len(reg_meta['source_files']),qc_policy=QC_POLICY,minimum_species=5,eligible_species=len(species),candidate_genes=len(genes),cached_methods=[m for m in METHODS if m!='tss100kb']),indent=2))
    print('Cached inputs ready:',len(genes),'genes',len(species),'species',flush=True)
    with ProcessPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(quantify_100kb,s) for s in selected]):future.result()
    for j,s in enumerate(species):consume(OUT/'tss100kb_cache'/f"{s['species']}.tsv.gz",j)
    np.savez_compressed(OUT/'inputs/tss100kb.npz',**base,y=values['tss100kb'],length_kb=length)
    print('All inputs ready in',round(time.time()-started,1),'seconds',flush=True)

if __name__=='__main__':main()
