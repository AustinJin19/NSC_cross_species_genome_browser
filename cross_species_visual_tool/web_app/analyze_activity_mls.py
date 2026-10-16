"""Genome-wide Pearson tests of unlogged gene-body+2kb activity versus raw MLS."""
import csv, json, math, argparse
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, false_discovery_control
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,default=ROOT/'output/gene_activity_signac_regions/species_median_by_symbol.tsv')
parser.add_argument('--output',type=Path,default=ROOT/'output/gene_activity_body_mls')
args=parser.parse_args()
out=args.output;out.mkdir(parents=True,exist_ok=True)
source=args.source
life=json.loads((Path(__file__).parent/'data/species_mls.json').read_text())
species={'NMR':'Heterocephalus glaber','BMR':'Nannospalax galili','DMR':'Fukomys damarensis','Mouse':'Mus musculus','Rat':'Rattus norvegicus','Rabbit':'Sylvilagus floridanus','ASM':'Acomys cahirinus'}
mls={k:life[v]['years'] for k,v in species.items()}
results=[];skipped=0
with source.open() as f:
 for row in csv.DictReader(f,delimiter='\t'):
  available=[s for s in species if row[s]!='NA' and math.isfinite(float(row[s]))]
  if len(available)<3:
   skipped+=1;continue
  x=np.array([mls[s] for s in available]);y=np.array([float(row[s]) for s in available])
  if np.ptp(x)==0 or np.ptp(y)==0:
   skipped+=1;continue
  r,p=pearsonr(x,y)
  slope=float(np.dot(x-x.mean(),y-y.mean())/np.dot(x-x.mean(),x-x.mean()))
  results.append(dict(gene=row['symbol_key'],n_species=len(available),pearson_r=float(r),p_value=float(p),slope=slope,species=','.join(available),**{s:row[s] for s in species}))
q=false_discovery_control([r['p_value'] for r in results],method='bh')
for r,v in zip(results,q):r['BH_FDR']=float(v)
results.sort(key=lambda r:r['p_value'])
fields=['gene','n_species','pearson_r','p_value','BH_FDR','slope','species',*species]
for name,rows in [('all_tested.tsv',results),('significant_FDR05.tsv',[r for r in results if r['BH_FDR']<.05]),('nominal_p05.tsv',[r for r in results if r['p_value']<.05]),('nominal_p05_positive.tsv',[r for r in results if r['p_value']<.05 and r['pearson_r']>0]),('nominal_p05_negative.tsv',[r for r in results if r['p_value']<.05 and r['pearson_r']<0])]:
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
summary=dict(source=str(source),nominal_positive=sum(r['p_value']<.05 and r['pearson_r']>0 for r in results),nominal_negative=sum(r['p_value']<.05 and r['pearson_r']<0 for r in results),tested=len(results),skipped=skipped,nominal_p05=sum(r['p_value']<.05 for r in results),FDR05=sum(r['BH_FDR']<.05 for r in results),MLS_years=mls,by_species_count={n:dict(tested=sum(r['n_species']==n for r in results),FDR05=sum(r['n_species']==n and r['BH_FDR']<.05 for r in results)) for n in range(3,8)},top=results[:10],top_complete=[r for r in results if r['n_species']==7][:10])
(out/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
