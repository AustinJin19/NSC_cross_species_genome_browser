"""Coverage and leave-one-species-out checks for expanded gene-body results."""
import csv,json,math
from pathlib import Path
import numpy as np
from scipy.stats import t
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/gene_activity_body_mls_expanded'
provenance=json.loads((OUT/'species_provenance.json').read_text());mls={r['species']:r['MLS'] for r in provenance};genes={}
for sid in mls:
 with (OUT/'species_scores'/f'{sid}.tsv').open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   r['activity']=float(r['activity']);r['length_kb']=float(r['length_kb'])
   old=genes.setdefault(r['gene'],{}).get(sid)
   if old is None or r['activity']>old['activity'] or (r['activity']==old['activity'] and r['locus']<old['locus']):genes[r['gene']][sid]=r
summary={}
for axis in ['MLS','log10_MLS']:
 rows=list(csv.DictReader((OUT/axis/'all_tested.tsv').open(),delimiter='\t'))
 for r in rows:
  for k in ['partial_r','p_value','BH_FDR']:r[k]=float(r[k])
  r['n_species']=int(r['n_species'])
 robust=[r for r in rows if r['n_species']>=20 and r['BH_FDR']<.05]
 with (OUT/axis/'FDR05_at_least_20_species.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(robust)
 selected=[r for r in rows if r['n_species']>=20][:30]
 checks=[]
 for r in selected:
  sp=genes[r['gene']];ids=sorted(sp);fits=[]
  for omit in [None]+ids:
   keep=[s for s in ids if s!=omit];x=np.array([mls[s] for s in keep]);x=np.log10(x) if axis=='log10_MLS' else x
   y=np.array([sp[s]['activity'] for s in keep]);z=np.array([sp[s]['length_kb'] for s in keep]);X=np.column_stack([np.ones(len(keep)),x,z]);beta=np.linalg.lstsq(X,y,rcond=None)[0]
   if np.linalg.matrix_rank(X)<3:continue
   residual=y-X@beta;df=len(y)-3;se=math.sqrt(max(0,float(residual@residual/df*np.linalg.inv(X.T@X)[1,1])))
   p=float(2*t.sf(abs(beta[1]/se),df)) if se else 0
   fits.append((omit,float(beta[1]),p))
   if omit is None:assert math.isclose(p,r['p_value'],rel_tol=1e-6,abs_tol=1e-12),(r['gene'],p,r['p_value'])
  loo=fits[1:];bow=next((v for v in loo if v[0]=='Balaena_mysticetus'),None)
  checks.append(dict(gene=r['gene'],n_species=r['n_species'],partial_r=r['partial_r'],BH_FDR=r['BH_FDR'],all_leave_one_out_same_sign=all((b>0)==(r['partial_r']>0) for _,b,p in loo),leave_one_out_max_p=max(p for _,b,p in loo),without_bowhead_p=bow[2] if bow else None))
 with (OUT/axis/'top30_coverage_sensitivity.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(checks[0]),delimiter='\t');w.writeheader();w.writerows(checks)
 summary[axis]={'at_least_20_species_tested':sum(r['n_species']>=20 for r in rows),'at_least_20_species_FDR05':len(robust),'at_least_20_species_FDR05_positive':sum(r['partial_r']>0 for r in robust),'at_least_20_species_FDR05_negative':sum(r['partial_r']<0 for r in robust),'top_positive':[r for r in robust if r['partial_r']>0][:8],'top_negative':[r for r in robust if r['partial_r']<0][:8],'sensitivity':checks[:8]}
(OUT/'coverage_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
