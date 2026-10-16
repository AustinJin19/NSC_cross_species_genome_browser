"""Peak-count regulatory strength versus lifespan, raw and log10 count sensitivity."""
import csv,json,sqlite3,statistics,math,argparse
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr,false_discovery_control
from activity import DB
from regulatory_strength import score_records
from server import SPECS
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/regulatory_strength';OUT.mkdir(parents=True,exist_ok=True)
life=json.loads((DB.parent/'species_mls.json').read_text());latin={s[0]:s[2] for s in SPECS}
con=sqlite3.connect(DB);con.row_factory=sqlite3.Row
tables={r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
records=[dict(r) for table in ['activity','activity_additions'] if table in tables for r in con.execute('SELECT * FROM '+table)]
for r in records:r['samples']=json.loads(r['samples'])
con.close()
def write(path,rows,fields):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
parser=argparse.ArgumentParser();parser.add_argument('--windows',nargs='+',type=int,choices=[20,50,100,200,500],default=[20,50,100,200,500]);args=parser.parse_args()
summary=json.loads((OUT/'summary.json').read_text()) if (OUT/'summary.json').exists() else {}
for flank in [w*1000 for w in args.windows]:
 folder=OUT/f'tss{flank//1000}kb';folder.mkdir(exist_ok=True);genes={};coverage=[]
 for sid in latin:
  rs=score_records([r for r in records if r['species']==sid],flank)
  write(folder/f'{sid}_locus_scores.tsv',[dict(gene=r['symbol_key'],locus=r['locus'],median=r['median'],samples=json.dumps(r['samples'])) for r in rs],['gene','locus','median','samples'])
  coverage.append(dict(species=sid,loci=len(rs),scored=sum(r['median'] is not None for r in rs)))
  if sid not in ['Human','Macaque','Balaena_mysticetus'] and life.get(latin[sid],{}).get('years'):
   for r in sorted(rs,key=lambda r:r['locus']):
    if r['median'] is None:continue
    old=genes.setdefault(r['symbol_key'],{}).get(sid)
    if old is None or r['median']>old:genes[r['symbol_key']][sid]=r['median']
  print(flank,sid,coverage[-1],flush=True)
 write(folder/'coverage.tsv',coverage,['species','loci','scored'])
 for scale in ['raw','log10']:
  for axis in ['MLS','log10_MLS']:
   results=[]
   for gene,values in genes.items():
    ids=sorted(s for s in values if scale=='raw' or values[s]>0)
    if len(ids)<5:continue
    x=np.array([life[latin[s]]['years'] for s in ids]);y=np.array([values[s] for s in ids]);x=np.log10(x) if axis=='log10_MLS' else x;y=np.log10(y) if scale=='log10' else y
    if np.ptp(x)==0 or np.ptp(y)==0:continue
    r,p=pearsonr(x,y);results.append(dict(gene=gene,n_species=len(ids),r=float(r),p_value=float(p),species=','.join(ids)))
   q=false_discovery_control([r['p_value'] for r in results],method='bh')
   for r,v in zip(results,q):r['BH_FDR']=float(v)
   results.sort(key=lambda r:r['p_value']);dest=folder/f'{scale}_{axis}';dest.mkdir(exist_ok=True)
   write(dest/'all_tested.tsv',results,list(results[0]))
   for sign in ['positive','negative']:
    sel=[r for r in results if r['BH_FDR']<.05 and (r['r']>0)==(sign=='positive')];write(dest/f'FDR05_{sign}.tsv',sel,list(results[0]))
   summary[f'{flank}_{scale}_{axis}']=dict(tested=len(results),positive=sum(r['BH_FDR']<.05 and r['r']>0 for r in results),negative=sum(r['BH_FDR']<.05 and r['r']<0 for r in results),top=results[:5])
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
