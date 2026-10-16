"""Genome-wide gene-body + 2kb MLS partial correlations controlling for length and body mass, matching browser rules."""
import csv,json,math,sqlite3,time,argparse
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import t, false_discovery_control
from server import SPECS,MODEL_INDEX,EXCLUDED_SAMPLES
from calculate_gene_activity import gene_regions,quantify
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--log-activity',action='store_true',help='Use log10 of the species median; exclude nonpositive values without a pseudocount')
parser.add_argument('--log-length',action='store_true',help='Adjust for log10 gene-body length rather than linear length')
parser.add_argument('--exclude-bowhead',action='store_true',help='Exclude Balaena mysticetus from this analysis')
parser.add_argument('--no-mass',action='store_true',help='Adjust for gene length only, without body mass')
args=parser.parse_args()
out_name='gene_activity_body_log10_mls_length_mass' if args.log_activity else 'gene_activity_body_mls_length_mass'
if args.log_length:out_name=out_name.replace('_length_mass','_loglength_mass')
if args.no_mass:out_name=out_name.replace('_mass','_no_mass')
if args.exclude_bowhead:out_name += '_no_bowhead'
OUT=ROOT/'output'/out_name;OUT.mkdir(exist_ok=True,parents=True)
CACHE=ROOT/'output/gene_activity_body_mls_expanded/species_scores';CACHE.mkdir(exist_ok=True)
LIFE=json.loads((Path(__file__).parent/'data/species_mls.json').read_text())
SELECTED=[s for s in SPECS if s[0] not in ('Human','Macaque') and LIFE.get(s[2],{}).get('years') and (args.no_mass or LIFE[s[2]].get('body_mass_g'))]
if args.exclude_bowhead:SELECTED=[s for s in SELECTED if s[0]!='Balaena_mysticetus']
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
  ids=sorted(s for s in sp if (not args.log_activity or sp[s]['activity']>0) and (not args.log_length or sp[s]['length_kb']>0));n=len(ids)
  if n<5:continue
  x=np.array([mls[s] for s in ids]);x=np.log10(x) if axis=='log10_MLS' else x
  y=np.array([sp[s]['activity'] for s in ids]);z=np.array([sp[s]['length_kb'] for s in ids])
  if args.log_activity:y=np.log10(y)
  if args.log_length:z=np.log10(z)
  mass=None if args.no_mass else np.array([math.log10(LIFE[next(v[2] for v in SELECTED if v[0]==sid)]['body_mass_g']) for sid in ids])
  C=np.column_stack([np.ones(n),z] if args.no_mass else [np.ones(n),z,mass])
  if np.linalg.matrix_rank(C)<C.shape[1]:continue
  rx=x-C@np.linalg.lstsq(C,x,rcond=None)[0];ry=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
  if rx@rx<=1e-12*np.sum((x-x.mean())**2) or ry@ry<=1e-12*np.sum((y-y.mean())**2):continue
  r=float(np.clip((rx@ry)/math.sqrt((rx@rx)*(ry@ry)),-1,1));df=n-C.shape[1]-1
  p=float(2*t.sf(abs(r)*math.sqrt(df/max(1-r*r,1e-300)),df))
  coefs=np.linalg.lstsq(np.column_stack([np.ones(n),x,z] if args.no_mass else [np.ones(n),x,z,mass]),y,rcond=None)[0]
  beta,bz=map(float,coefs[1:3])
  record=dict(gene=gene,n_species=n,partial_r=r,p_value=p,MLS_coefficient=beta)
  record['length_coefficient_per_log10kb' if args.log_length else 'length_coefficient_per_kb']=bz
  if not args.no_mass:record['mass_coefficient_per_log10g']=float(coefs[3])
  record.update(df=df,species=','.join(ids))
  results.append(record)
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
