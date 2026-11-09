"""Compare gene-body accessibility and liver RNA correlations against raw MLS."""
import csv, json, sqlite3, math
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import pearsonr, false_discovery_control

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
DATA=ROOT/'cross_species_visual_tool/web_app/data'
species={'NMR':'Heterocephalus glaber','BMR':'Nannospalax galili','DMR':'Fukomys damarensis','Mouse':'Mus musculus','Rat':'Rattus norvegicus','Rabbit':'Sylvilagus floridanus','ASM':'Acomys cahirinus'}
reverse={v:k for k,v in species.items()}
life=json.loads((DATA/'species_mls.json').read_text())
mls={k:life[v]['years'] for k,v in species.items()}
body={r['gene']:r for r in csv.DictReader((ROOT/'output/gene_activity_body_mls/all_tested.tsv').open(),delimiter='\t')}
rna=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
con=sqlite3.connect(f'file:{DATA}/liver_omics.sqlite?mode=ro',uri=True)
for gene,sp,raw,log in con.execute("SELECT symbol,species,raw,log FROM measurements WHERE assay='rna' AND sample!='BMR_Turk_Lu'"):
 if sp not in reverse or gene not in body: continue
 for metric,value in [('log',log),('raw',raw)]:
  if value is not None and math.isfinite(value):rna[metric][gene][reverse[sp]].append(value)
con.close()

def counts(rows):
 pos=sum(r['atac_r']>0 and r['rna_r']>0 for r in rows)
 neg=sum(r['atac_r']<0 and r['rna_r']<0 for r in rows)
 return dict(total=len(rows),positive_both=pos,negative_both=neg,same=pos+neg,opposite=sum(r['atac_r']*r['rna_r']<0 for r in rows),same_percent=100*(pos+neg)/len(rows) if rows else None)

summaries={}
for metric,mode in [('log','matched'),('raw','matched'),('log','original_body'),('log','matched_without_ASM')]:
 rows=[]
 for gene,b in body.items():
  avail=[s for s in species if rna[metric][gene].get(s) and b[s]!='NA' and math.isfinite(float(b[s])) and not(mode=='matched_without_ASM' and s=='ASM')]
  if len(avail)<3:continue
  x=np.array([mls[s] for s in avail]);a=np.array([float(b[s]) for s in avail]);y=np.array([np.median(rna[metric][gene][s]) for s in avail])
  if np.ptp(a)==0 or np.ptp(y)==0:continue
  ar,ap=pearsonr(x,a);rr,rp=pearsonr(x,y)
  if mode=='original_body':ar,ap=float(b['pearson_r']),float(b['p_value'])
  assert math.isclose(rr,float(np.corrcoef(x,y)[0,1]),abs_tol=1e-12)
  rows.append(dict(gene=gene,n_common=len(avail),species=','.join(avail),atac_r=float(ar),atac_p=float(ap),rna_r=float(rr),rna_p=float(rp),direction='positive_both' if ar>0 and rr>0 else 'negative_both' if ar<0 and rr<0 else 'opposite',original_atac_r=float(b['pearson_r']),original_atac_p=float(b['p_value']),original_atac_FDR=float(b['BH_FDR']),**{f'{s}_ATAC':float(b[s]) if s in avail else 'NA' for s in species},**{f'{s}_RNA':float(np.median(rna[metric][gene][s])) if s in avail else 'NA' for s in species}))
 for assay in ['atac','rna']:
  for r,q in zip(rows,false_discovery_control([r[assay+'_p'] for r in rows])):r[assay+'_FDR']=float(q)
 groups={'all':rows,'atac_p05':[r for r in rows if r['atac_p']<.05],'rna_p05':[r for r in rows if r['rna_p']<.05],'both_p05':[r for r in rows if r['atac_p']<.05 and r['rna_p']<.05],'both_FDR05':[r for r in rows if r['atac_FDR']<.05 and r['rna_FDR']<.05]}
 key=mode+'_'+metric
 summaries[key]={k:counts(v) for k,v in groups.items()}
 summaries[key]['by_n']={n:counts([r for r in rows if r['n_common']==n]) for n in range(3,7)}
 for name,subset in [('all',rows),('same_direction',[r for r in rows if r['atac_r']*r['rna_r']>0]),('same_direction_both_p05',[r for r in groups['both_p05'] if r['atac_r']*r['rna_r']>0])]:
  subset.sort(key=lambda r:max(r['atac_p'],r['rna_p']))
  with (OUT/f'{key}_{name}.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(subset)
(OUT/'summary.json').write_text(json.dumps(summaries,indent=2))
(OUT/'METHODS.txt').write_text('Gene-body + strand-aware 2kb promoter unlogged bigWig sum, species medians, versus liver RNA species medians. Pearson correlations against raw MLS years. Primary RNA metric: source log.tpm; raw TPM sensitivity included. RNA from existing liver_omics.sqlite (CrossSpecies_Count primary, rnaLongTPM mouse supplement). BMR_Turk_Lu excluded. Human and macaque excluded. Rabbit absent from RNA. Primary analysis uses the identical available species in both assays for each gene, minimum three, nonconstant values. Original-body comparison retains original ATAC coefficients (potentially including Rabbit) while RNA uses overlapping species. BH correction separately within each assay across the comparable genes; original genome-wide ATAC FDR retained separately. Symbols matched case-insensitively via existing indexes, not verified orthologs. No phylogenetic adjustment. ASM mapping uses source species field Acomys cahirinus; some source common labels say Acomys russatus, so an ASM-excluded sensitivity is included. No source data or website modified.\n')
print(json.dumps(summaries,indent=2))
