"""Check the live site's sample exclusions and agreement with batch outputs."""
import concurrent.futures, csv, gzip, json, math, statistics, sys, urllib.request
from pathlib import Path
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT/'cross_species_visual_tool/web_app'))
from sample_qc import EXCLUDED_SAMPLES,sample_id

def get(path):
    with urllib.request.urlopen('http://127.0.0.1:8765'+path,timeout=180) as r:return json.load(r)
catalog=get('/api/catalog')
assert len(catalog)==33 and sum(len(s['samples']) for s in catalog)==77
assert not any(sample_id(x) in EXCLUDED_SAMPLES for s in catalog for x in s['samples'])
assert not any(s['id']=='Balaena_mysticetus' for s in catalog)
assert next(s for s in catalog if s['id']=='ASM')['samples']==['ASM3Liver','ASM4Liver']
strategies=['body','tss2kb','tss5kb','tss100kb','tss200kb','tss500kb']+[f'{prefix}{k}kb' for prefix in ['peaks','peak_signal'] for k in [20,50,100,200,500]]
jobs=[('SIRT6',s) for s in strategies]+[(g,s) for g in ['IGF1','HAS2'] for s in ['body','tss2kb','tss5kb','tss200kb','tss500kb','peaks20kb','peaks50kb','peaks100kb','peaks200kb','peaks500kb']]
def check(job):
    gene,strategy=job; data=get('/api/gene-activity?gene='+gene+'&strategy='+strategy)
    chosen={}
    for r in data['records']:
        assert r['species']!='Balaena_mysticetus'
        assert not any(sample_id(x) in EXCLUDED_SAMPLES for x in r['samples'])
        assert not any(sample_id(x) in EXCLUDED_SAMPLES for x in r.get('window_details',{}))
        values=[v for v in r['samples'].values() if v is not None and math.isfinite(v)]
        assert r['median']==(statistics.median(values) if values else None)
        if r['median'] is None:continue
        old=chosen.get(r['species'])
        if old is None or (-r['median'],r['locus'])<(-old['median'],old['locus']):chosen[r['species']]=r
    return job,chosen
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    live=dict(pool.map(check,jobs))
    print('Live API: all 36 queries passed sample and median checks.',flush=True)
comparisons=0
max_relative_signal_difference=0.0
for folder in ['activity_selected_loci','regulatory_selected_loci']:
    for path in (OUT/folder).glob('*.tsv.gz'):
        with gzip.open(path,'rt') as f:
            for row in csv.DictReader(f,delimiter='\t'):
                key=(row['gene'],row['method'])
                if key not in live:continue
                actual=live[key].get(row['species'])
                assert actual is not None,(key,row['species'])
                assert actual['locus']==row['locus'],(key,row['species'],actual['locus'],row['locus'])
                assert math.isclose(actual['median'],float(row['median']),rel_tol=1e-6,abs_tol=1e-6),(key,row['species'],actual['median'],row['median'])
                assert actual['region_eligible']==(row['region_eligible']=='True')
                assert math.isclose(actual['available_fraction'],float(row['available_fraction']),abs_tol=1e-12)
                # Existing per-sample signal caches can have float32 rounding.
                max_relative_signal_difference=max(max_relative_signal_difference,abs(actual['median']-float(row['median']))/max(abs(actual['median']),1e-12))
                comparisons+=1
report=dict(species=33,samples=77,excluded_samples=len(EXCLUDED_SAMPLES),api_queries=len(jobs),batch_locus_comparisons=comparisons,max_relative_signal_difference=max_relative_signal_difference,relative_signal_tolerance=1e-6,passed=True)
(OUT/'live_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
