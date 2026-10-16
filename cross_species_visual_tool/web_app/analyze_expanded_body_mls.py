"""Genome-wide gene-body + 2kb MLS partial correlations, matching browser rules."""
import csv,json,math,sqlite3,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import t, false_discovery_control
from server import SPECS,MODEL_INDEX,EXCLUDED_SAMPLES
from calculate_gene_activity import gene_regions,quantify
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/gene_activity_body_mls_expanded';OUT.mkdir(exist_ok=True,parents=True)
CACHE=OUT/'species_scores';CACHE.mkdir(exist_ok=True)
LIFE=json.loads((Path(__file__).parent/'data/species_mls.json').read_text())
SELECTED=[s for s in SPECS if s[0] not in ('Human','Macaque') and LIFE.get(s[2],{}).get('years')]
def write(path,rows,fields):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def length(locus):
 a,b=map(int,locus.rsplit('|',1)[1].rsplit(':',2)[1].split('-'));return (b-a)/1000
con=sqlite3.connect(Path(__file__).parent/'data/gene_activity.sqlite');con.row_factory=sqlite3.Row
existing={r[0] for r in con.execute('SELECT DISTINCT species FROM activity')}
genes=defaultdict(dict);provenance=[]
for sid,name,latin,color,folder,prefix in SELECTED:
 start=time.time();cache=CACHE/(sid+'.tsv')
 if cache.exists():
  rows=list(csv.DictReader(cache.open(),delimiter='\t'))
  for r in rows:r['activity']=float(r['activity']);r['length_kb']=float(r['length_kb'])
 elif sid in existing:
  rows=[dict(gene=r['symbol_key'],locus=r['locus'],activity=r['median'],length_kb=length(r['locus'])) for r in con.execute('SELECT * FROM activity WHERE species=?',(sid,)) if r['median'] is not None and math.isfinite(r['median'])]
  write(cache,rows,['gene','locus','activity','length_kb'])
 else:
  regions=gene_regions(ROOT/'cross_species_visual_tool/index'/f'{MODEL_INDEX.get(sid,sid)}_models.tsv.gz')
  paths=sorted(p for p in (ROOT/'cross_species_visual_tool'/folder).glob('*.bw') if p.name.lower().startswith(prefix.lower()) and p.stem.removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES)
  arrays=[]
  for path in paths:
   values,*_=quantify(path,regions,chunked=False);arrays.append(values)
   print(sid,path.name,'quantified',flush=True)
  matrix=np.asarray(arrays).T
  rows=[]
  for region,values in zip(regions,matrix):
   finite=values[np.isfinite(values)]
   if len(finite):rows.append(dict(gene=region['symbol_key'],locus=region['locus_id'],activity=float(np.median(finite)),length_kb=length(region['locus_id'])))
  write(cache,rows,['gene','locus','activity','length_kb'])
 for r in sorted(rows,key=lambda r:r['locus']):
  old=genes[r['gene']].get(sid)
  if old is None or r['activity']>old['activity']:genes[r['gene']][sid]=r
 provenance.append(dict(species=sid,scientific_name=latin,MLS=LIFE[latin]['years'],MLS_source=LIFE[latin],loci=len(rows)))
 print(sid,len(rows),'loci',round(time.time()-start,1),'seconds',flush=True)
con.close()
mls={s[0]:LIFE[s[2]]['years'] for s in SELECTED}
summary={}
for axis in ['MLS','log10_MLS']:
 results=[]
 for gene,sp in genes.items():
  ids=sorted(sp);n=len(ids)
  if n<4:continue
  x=np.array([mls[s] for s in ids]);x=np.log10(x) if axis=='log10_MLS' else x
  y=np.array([sp[s]['activity'] for s in ids]);z=np.array([sp[s]['length_kb'] for s in ids])
  xc=x-x.mean();yc=y-y.mean();zc=z-z.mean();zz=zc@zc;xx=xc@xc;yy=yc@yc
  if min(zz,xx,yy)<=0:continue
  rx=xc-(xc@zc)/zz*zc;ry=yc-(yc@zc)/zz*zc
  if rx@rx<=1e-12*xx or ry@ry<=1e-12*yy:continue
  r=float(np.clip((rx@ry)/math.sqrt((rx@rx)*(ry@ry)),-1,1));df=n-3
  p=float(2*t.sf(abs(r)*math.sqrt(df/max(1-r*r,1e-300)),df))
  beta=float((rx@ry)/(rx@rx));bz=float((yc@zc-beta*(xc@zc))/zz)
  results.append(dict(gene=gene,n_species=n,partial_r=r,p_value=p,MLS_coefficient=beta,length_coefficient_per_kb=bz,df=df,species=','.join(ids)))
 q=false_discovery_control([r['p_value'] for r in results],method='bh')
 for r,v in zip(results,q):r['BH_FDR']=float(v)
 results.sort(key=lambda r:r['p_value']);fields=list(results[0]);folder=OUT/axis;folder.mkdir(exist_ok=True)
 write(folder/'all_tested.tsv',results,fields)
 for threshold,col in [('FDR05','BH_FDR'),('nominal_p05','p_value')]:
  for direction in ['positive','negative']:
   rows=[r for r in results if r[col]<.05 and ((r['partial_r']>0)==(direction=='positive'))]
   write(folder/(threshold+'_'+direction+'.tsv'),rows,fields)
 summary[axis]={'tested':len(results),'nominal_positive':sum(r['p_value']<.05 and r['partial_r']>0 for r in results),'nominal_negative':sum(r['p_value']<.05 and r['partial_r']<0 for r in results),'FDR_positive':sum(r['BH_FDR']<.05 and r['partial_r']>0 for r in results),'FDR_negative':sum(r['BH_FDR']<.05 and r['partial_r']<0 for r in results),'top_positive':[r for r in results if r['partial_r']>0][:10],'top_negative':[r for r in results if r['partial_r']<0][:10]}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(OUT/'species_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
