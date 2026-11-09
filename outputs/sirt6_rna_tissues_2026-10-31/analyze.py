"""Tissue-specific RNA rank scores and pooled-lambda PGLS for the SIRT6 locus."""
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import false_discovery_control, rankdata, t

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
(OUT/'build').mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT/'outputs/pgls_genomewide_2026-10-27'))
from pgls_batch import fit_batch
prior = json.loads((ROOT/'outputs/sirt6_shared_lambda_rank_scores_2026-10-31/analysis.json').read_text())
matrix = json.loads((ROOT/'outputs/sirt6_neighborhood_rank_matrices_2026-10-29/matrix_data.json').read_text())
genes = matrix['metadata']['genes']
species = prior['species']
included_species = [s for s in species if s['group'] != 'Excluded']
by_latin = {s['latin']:s for s in species}
aliases = {'Ellobius lutescens':'Bramus lutescens'}
sources = [json.loads((OUT/(name+'_extracted.json')).read_text()) for name in ['CrossSpecies_Count','rnaLongTPM']]

# Match the established browser source priority: primary workbook, then the
# mouse supplement. Audit any genuinely new non-mouse records before exclusion.
chosen = {}
duplicate_checks = []
new_nonmouse = []
for idx, source in enumerate(sources):
    local_seen = set()
    for r in source['records']:
        latin = aliases.get(r['species'], r['species'])
        if latin not in by_latin: continue
        key = (latin, r['tissue'], r['sample'], r['gene'])
        assert key not in local_seen, ('Duplicate source measurement', key, source['workbook'])
        local_seen.add(key)
        if key in chosen:
            old = chosen[key]
            same = all((old[k] is None and r[k] is None) or
                (old[k] is not None and r[k] is not None and np.isclose(old[k],r[k],rtol=1e-10,atol=1e-12))
                for k in ['tpm','log_tpm'])
            duplicate_checks.append(dict(tissue=r['tissue'],gene=r['gene'],sample=r['sample'],
                species=latin,values_agree=bool(same),retained=old['workbook'],other=r['workbook']))
            continue
        if idx and latin != 'Mus musculus':
            new_nonmouse.append(r)
            continue
        chosen[key] = dict(r, latin=latin, species_id=by_latin[latin]['id'])
assert not new_nonmouse, f'{len(new_nonmouse)} new non-mouse rows require a source-selection review'
records = list(chosen.values())
groups = defaultdict(list)
for r in records:
    if r['log_tpm'] is not None:
        groups[(r['tissue'],r['species_id'],r['gene'])].append(r['log_tpm'])
values = {key:float(np.median(v)) for key,v in groups.items()}
tissues = sorted({r['tissue'] for r in records}, key=lambda x:(x!='Liver',x))

# The liver re-extraction must reproduce the previous analysis without changes.
old_values = {(r['species_id'],r['gene']):r['value'] for r in prior['measurements'] if r['assay']=='rna'}
for key, old in old_values.items():
    got = values.get(('Liver',*key))
    assert (got is None and old is None) or (got is not None and old is not None and abs(got-old)<1e-12), (key,got,old)

tree = json.loads((ROOT/'cross_species_visual_tool/web_app/data/phylogenetic_covariance.json').read_text())
index = {s:i for i,s in enumerate(tree['species'])}
full_C = np.asarray(tree['covariance'])
results, measurements, prepared, profiles = [], [], [], []
panels = {}
for tissue in tissues:
    for gene in genes:
        eligible = sorted([s for s in included_species if (tissue,s['id'],gene) in values], key=lambda s:s['latin'])
        common = [g for g in genes if eligible and all((tissue,s['id'],g) in values for s in eligible)]
        panels[tissue+'|'+gene] = common
        scores = {}
        for s in eligible:
            rank = float(rankdata([-values[tissue,s['id'],g] for g in common],method='average')[common.index(gene)])
            scores[s['id']] = (rank, 1-(rank-1)/(len(common)-1) if len(common)>1 else None)
        counts = Counter(s['group'] for s in eligible)
        r = dict(tissue=tissue,gene=gene,status='Tested',n=len(eligible),
            n_rearranged=counts['Rearranged'],n_unrearranged=counts['Unrearranged'],reference_genes=len(common),
            lambda_ml=None,beta=None,se=None,ci_low=None,ci_high=None,p=None,q_tissue=None,q_global=None,
            q_sirt6_tissues=None,intercept=None,df=None,log_likelihood=None,gls_r2=None,unique_scores=None,
            median_rearranged=None,median_unrearranged=None,bm_beta=None,bm_p=None)
        if min(r['n_rearranged'],r['n_unrearranged'])<3:
            r['status']='Not tested: fewer than 3 species in a group'
        elif len(common)<2:
            r['status']='Not tested: fewer than 2 shared ranking genes'
        else:
            y = np.asarray([scores[s['id']][1] for s in eligible])
            x = np.asarray([int(s['group']=='Rearranged') for s in eligible])
            r.update(unique_scores=len(np.unique(y)),median_rearranged=float(np.median(y[x==1])),
                median_unrearranged=float(np.median(y[x==0])))
            if np.ptp(y)<1e-14:
                r['status']='Not estimable: constant rank score'
            else:
                ids = [index[s['latin']] for s in eligible]
                C = full_C[np.ix_(ids,ids)]
                X = np.column_stack([np.ones(len(y)),x])
                f = fit_batch(y[None,:],X[None,:,:],C,0)
                if not f['valid'][0]:
                    r['status']='Not estimable: zero residual variance'
                else:
                    prepared.append(dict(result=r,Y=y[None,:],X=X[None,:,:],C=C,species=eligible))
        results.append(r)
        for s in species:
            rank, score = scores.get(s['id'],(None,None))
            measurements.append(dict(tissue=tissue,gene=gene,species_id=s['id'],species=s['common'],latin=s['latin'],
                group=s['group'],median_log_tpm=values.get((tissue,s['id'],gene)),
                n_samples=len(groups.get((tissue,s['id'],gene),[])),rank=rank,reference_genes=len(common),
                rank_score=score,rearranged=int(s['group']=='Rearranged') if s['group']!='Excluded' else None,
                used_in_test=r['status']=='Tested' and s['id'] in scores))

def fit_one(p,lam):
    f = fit_batch(p['Y'],p['X'],p['C'],lam)
    assert f['valid'][0]
    return f

pools = {}
cases = []
for tissue in tissues:
    entries = [p for p in prepared if p['result']['tissue']==tissue]
    if not entries:
        pools[tissue] = dict(lambda_value=None,genes=0,sum_profile_log_likelihood=None)
        continue
    objective = lambda lam:sum(float(fit_one(p,lam)['log_likelihood'][0]) for p in entries)
    grid = np.linspace(0,1,101)
    ll = np.asarray([objective(lam) for lam in grid])
    j = int(ll.argmax()); lam = float(grid[j]); best = float(ll[j])
    for i in range(len(grid)):
        if (i and ll[i]<ll[i-1]) or (i<len(grid)-1 and ll[i]<ll[i+1]):continue
        z = minimize_scalar(lambda x:-objective(x),bounds=(max(0,grid[i]-.01),min(1,grid[i]+.01)),
            method='bounded',options={'xatol':1e-11})
        if -z.fun>best:lam,best=float(z.x),float(-z.fun)
    pools[tissue] = dict(lambda_value=lam,genes=len(entries),sum_profile_log_likelihood=best)
    profiles.extend(dict(tissue=tissue,lambda_value=float(l),log_likelihood=float(v),relative_to_max=float(v-best)) for l,v in zip(grid,ll))
    for p in entries:
        f = fit_one(p,lam); r = p['result']; bm = fit_one(p,1)
        beta,se,df = float(f['beta'][0,0]),float(f['se'][0,0]),int(f['df'][0])
        crit = t.ppf(.975,df)
        r.update(lambda_ml=lam,beta=beta,se=se,df=df,ci_low=float(beta-crit*se),ci_high=float(beta+crit*se),
            p=float(f['p'][0]),intercept=float(f['intercept'][0]),log_likelihood=float(f['log_likelihood'][0]),
            gls_r2=float(f['r_squared'][0]),bm_beta=float(bm['beta'][0,0]),bm_p=float(bm['p'][0]))
        # Direct GLS is independent of the spectral-whitening implementation.
        V=p['C']*lam;np.fill_diagonal(V,np.diag(p['C']))
        X=p['X'][0]; y=p['Y'][0]
        viX=np.linalg.solve(V,X); viy=np.linalg.solve(V,y)
        coef=np.linalg.solve(X.T@viX,X.T@viy);resid=y-X@coef
        se_check=np.sqrt(np.diag(np.linalg.inv(X.T@viX))*(resid@np.linalg.solve(V,resid))/df)
        assert abs(coef[1]-beta)<1e-10 and abs(se_check[1]-se)<1e-10
        assert abs(2*t.sf(abs(coef[1]/se_check[1]),df)-r['p'])<1e-10
        cases.append(dict(assay=tissue,gene=r['gene'],points=[dict(species=s['latin'].replace(' ','_'),
            x=int(s['group']=='Rearranged'),y=float(p['Y'][0,i])) for i,s in enumerate(p['species'])],
            covariance=p['C'].tolist(),expected=r))
    print(tissue, pools[tissue], flush=True)

def adjust(rows,key):
    rows=[r for r in rows if r['p'] is not None]
    if rows:
        for r,q in zip(rows,false_discovery_control([r['p'] for r in rows],method='bh')):r[key]=float(q)
for tissue in tissues:adjust([r for r in results if r['tissue']==tissue],'q_tissue')
adjust(results,'q_global')
adjust([r for r in results if r['gene']=='SIRT6'],'q_sirt6_tissues')

for old in prior['results']:
    if old['assay']!='rna':continue
    new=next(r for r in results if r['tissue']=='Liver' and r['gene']==old['gene'])
    assert new['reference_genes']==old['reference_genes']
    if old['p'] is not None:
        for key,old_key in [('p','p'),('beta','beta'),('q_tissue','q'),('lambda_ml','lambda_ml')]:
            assert abs(new[key]-old[old_key])<1e-10,(old['gene'],key,new[key],old[old_key])

summary=[]
for tissue in tissues:
    rr=[r for r in results if r['tissue']==tissue]
    available={r['species_id'] for r in records if r['tissue']==tissue and by_latin[r['latin']]['group']!='Excluded'}
    samples={(r['species_id'],r['sample']) for r in records if r['tissue']==tissue and r['species_id'] in available}
    summary.append(dict(tissue=tissue,species=len(available),samples=len(samples),
        rearranged=sum(s['id'] in available and s['group']=='Rearranged' for s in species),
        unrearranged=sum(s['id'] in available and s['group']=='Unrearranged' for s in species),
        tested=sum(r['p'] is not None for r in rr),constant=sum(r['status']=='Not estimable: constant rank score' for r in rr),
        nominal=sum(r['p'] is not None and r['p']<.05 for r in rr),
        fdr_tissue=sum(r['q_tissue'] is not None and r['q_tissue']<.05 for r in rr),
        fdr_global=sum(r['q_global'] is not None and r['q_global']<.05 for r in rr),lambda_value=pools[tissue]['lambda_value']))

sample_rows=[]
for s in sources:
    for r in s['samples']:
        latin=aliases.get(r['species'],r['species']);sp=by_latin.get(latin)
        if sp is None:continue
        selected=[z for z in records if z['workbook']==s['workbook'] and z['sample']==r['sample'] and z['latin']==latin and z['tissue']==r['tissue']]
        sample_rows.append(dict(workbook=s['workbook'],sample=r['sample'],tissue=r['tissue'],source_tissue=r['source_tissue'],
            species=sp['common'],latin=latin,group=sp['group'],source=r['source'],total_source_rows=r['rows'],
            selected_locus_genes=len(selected),used=bool(selected) and sp['group']!='Excluded'))

notes=[
    'Question: association between the SIRT6-locus rearrangement group and within-species RNA expression rank in each tissue. The predictor is rearrangement, not MLS.',
    'The same original 63 locus labels and 10 rearranged / 20 unrearranged species are retained. Human, macaque, bowhead whale and muskrat are excluded. User-resolved group assignments are unchanged.',
    'CrossSpecies_Count.xlsx supplies primary RNA measurements; rnaLongTPM.xlsx supplies missing mouse records. Duplicate sample/gene/tissue/species measurements are not counted twice. Original source records remain unchanged.',
    'BMR_Turk_Lu is assigned to Lung, following the previously confirmed correction; both input workbooks label it Liver. Source tissue and corrected tissue are retained.',
    'Use the species median of source log.tpm within tissue, without another logarithm. Each species contributes one observation, irrespective of sample number. Symbols are matched case-insensitively, with no new alias or orthology substitutions.',
    'For each tissue and target gene, use classified species with a measurement and rank only candidate genes measured in all of those species. Rank 1 is highest; ties receive average ranks. Rank score = 1 - (rank-1)/(K-1). Missing values are never zero.',
    'Reference gene sets and species coverage can differ by target and tissue. Therefore a difference in statistical significance between tissues is not itself a tissue-by-rearrangement interaction test.',
    'At least three measured species in each group and at least two shared ranking genes are required. Constant outcomes or zero residual variance have no estimable Gaussian test.',
    'Model: rank score ~ rearranged. One Pagel lambda is shared across all valid genes within each tissue, maximizing their summed full-model Gaussian profile log likelihood over [0,1]. Each gene has its own intercept, coefficient and residual variance. Grid step 0.01 plus local refinement and both endpoints.',
    'The pooled criterion is a composite likelihood because genes and their ranks are dependent. Conventional GLS t intervals and p-values use n-2 degrees of freedom conditional on the estimated lambda; lambda uncertainty is not propagated.',
    'TimeTree covariance uses the original full-tree root. Diagonal entries are root-to-tip lengths; off-diagonal entries are lambda times shared root-to-MRCA lengths. Lambda=0 gives ordinary regression here; lambda=1 is also supplied as a sensitivity check.',
    'BH q within tissue adjusts all valid gene tests in that tissue. BH q across tissues adjusts all valid gene-by-tissue tests jointly. SIRT6-only q adjusts SIRT6 tests across tested tissues and is clearly labeled separately.',
    'All rearranged species form one clade among the included taxa; the comparison does not establish inversion causation or independent evolutionary replication. Bounded, tied rank scores make Gaussian inference exploratory, especially when K is small.',
    'African spiny mouse source species is Acomys cahirinus while its source common-name field sometimes says Acomys russatus. The established species assignment is retained.',
    'Liver medians, shared panels and valid primary fits reproduce the previous liver-only shared-lambda analysis exactly. Website data and prior analysis outputs are unchanged.',
]
payload=dict(analysis_date='2026-10-31',genes=genes,tissues=tissues,species=species,results=results,summary=summary,
    measurements=measurements,panels=panels,shared_lambda=pools,profiles=profiles,samples=sample_rows,notes=notes,
    sources=[{k:v for k,v in s.items() if k not in ['samples','records']} for s in sources],
    duplicate_audit=dict(overlapping_records=len(duplicate_checks),discordant=sum(not r['values_agree'] for r in duplicate_checks)),
    tree_file=prior['tree_file'],tree_matrix=prior['tree_matrix'],tree_sha256=prior['tree_sha256'],references=prior['references'])
(OUT/'analysis.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
(OUT/'validation_cases.json').write_text(json.dumps(cases,indent=2,allow_nan=False)+'\n')
def tsv(name,rows):
    if not rows:return
    with (OUT/name).open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader()
        writer.writerows({k:('NA' if v is None else v) for k,v in row.items()} for row in rows)
tsv('RNA_tissue_gene_tests.tsv',results)
tsv('RNA_tissue_summary.tsv',summary)
tsv('RNA_tissue_rank_scores.tsv',measurements)
tsv('RNA_tissue_source_samples.tsv',sample_rows)
tsv('RNA_tissue_source_records.tsv',records)
tsv('RNA_tissue_duplicate_audit.tsv',duplicate_checks)
print(json.dumps(dict(summary=summary,SIRT6=[r for r in results if r['gene']=='SIRT6'],duplicates=payload['duplicate_audit']),indent=2),flush=True)

os.environ.setdefault('MPLCONFIGDIR',str(OUT/'build/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
valid=[r for r in results if r['gene']=='SIRT6']
if valid:
    fig,ax=plt.subplots(figsize=(8,max(3.5,len(valid)*.7+1.5)),layout='constrained')
    for i,r in enumerate(valid):
        if r['p'] is not None:
            color='#0072B2' if r['q_global']<.05 else '#82919E'
            ax.errorbar(r['beta'],i,xerr=[[r['beta']-r['ci_low']],[r['ci_high']-r['beta']]],fmt='o',c=color,capsize=4,ms=7)
            label=f"p={r['p']:.3g}; q={r['q_global']:.3g}; λ={r['lambda_ml']:.3f}"
        else:
            label='Constant rank: 2 of 3 in every species' if r['unique_scores']==1 else 'Not tested: only one species'
        ax.text(1.02,i,label,transform=ax.get_yaxis_transform(),va='center',fontsize=9)
    ax.set_yticks(range(len(valid)),[f"{r['tissue']} (n={r['n']}, K={r['reference_genes']})" for r in valid]);ax.invert_yaxis()
    ax.axvline(0,color='#333333',ls='--',lw=1);ax.set_xlabel('SIRT6 rank-score difference and 95% interval\nPositive = higher relative expression in rearranged species')
    ax.set_title('SIRT6 RNA expression ranks across tissues\nq adjusts all tested gene–tissue pairs',fontsize=12)
    ax.set_axisbelow(True);ax.xaxis.grid(color='#E5E7EB')
    for ext in ['png','svg','pdf']:fig.savefig(OUT/f'SIRT6_RNA_tissue_effects.{ext}',dpi=220,facecolor='white',bbox_inches='tight')
    plt.close(fig)
tested_genes=[g for g in genes if any(r['gene']==g and r['p'] is not None for r in results)]
tested_tissues=[s['tissue'] for s in summary if s['tested']]
lookup={(r['tissue'],r['gene']):r for r in results}
heat=np.full((len(tested_genes),len(tested_tissues)),np.nan)
for i,g in enumerate(tested_genes):
    for j,tissue in enumerate(tested_tissues):
        v=lookup[tissue,g]['beta']
        if v is not None:heat[i,j]=v
if heat.size:
    fig,ax=plt.subplots(figsize=(max(5,len(tested_tissues)*1.2),max(5,len(tested_genes)*.23)),layout='constrained')
    cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#E7E9EC')
    lim=max(.1,float(np.nanmax(abs(heat))))
    im=ax.imshow(heat,cmap=cm,vmin=-lim,vmax=lim,aspect='auto')
    ax.set_xticks(range(len(tested_tissues)),tested_tissues,rotation=35,ha='right');ax.set_yticks(range(len(tested_genes)),tested_genes,fontsize=8)
    for i,g in enumerate(tested_genes):
        for j,tissue in enumerate(tested_tissues):
            r=lookup[tissue,g]
            if r['q_global'] is not None and r['q_global']<.05:ax.text(j,i,'**',ha='center',va='center',fontsize=10)
            elif r['q_tissue'] is not None and r['q_tissue']<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=10)
    ax.set_title('RNA rank-score difference: rearranged minus unrearranged\nNo associations pass FDR < 0.05; gray = no estimable fit',fontsize=10)
    fig.colorbar(im,ax=ax,shrink=.5,label='Group coefficient')
    for ext in ['png','svg','pdf']:fig.savefig(OUT/f'RNA_locus_tissue_heatmap.{ext}',dpi=220,facecolor='white',bbox_inches='tight')
    plt.close(fig)
