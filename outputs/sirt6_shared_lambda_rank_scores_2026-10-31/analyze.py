"""Shared-lambda rank-score PGLS, estimated across genes separately by assay.

Summed univariate profile log likelihoods form a composite likelihood because
genes and within-species ranks are dependent. It is not a multivariate fit.
"""
import copy
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t, false_discovery_control

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
PREVIOUS=ROOT/'outputs/sirt6_pgls_rank_scores_2026-10-31'
sys.path.insert(0,str(ROOT/'outputs/pgls_genomewide_2026-10-27'))
from pgls_batch import fit_batch
(OUT/'build').mkdir(exist_ok=True)
prior=json.loads((PREVIOUS/'analysis.json').read_text())
cases=json.loads((PREVIOUS/'validation_cases.json').read_text())
payload=copy.deepcopy(prior)
results=payload['results']
result_lookup={(r['assay'],r['gene']):r for r in results}
prepared=[]
for c in cases:
    y=np.asarray([p['y'] for p in c['points']])
    X=np.column_stack([np.ones(len(y)),[p['x'] for p in c['points']]])
    C=np.asarray(c['covariance'])
    prepared.append(dict(case=c,Y=y[None,:],X=X[None,:,:],C=C))

def fit_case(p,lam):
    f=fit_batch(p['Y'],p['X'],p['C'],lambda_mode=lam)
    assert f['valid'][0]
    return f

def pooled_loglik(lam,entries):
    return sum(float(fit_case(p,lam)['log_likelihood'][0]) for p in entries)

pools={}
profiles=[]
for assay in ['rna','atac']:
    entries=[p for p in prepared if p['case']['assay']==assay]
    grid=np.linspace(0,1,101)
    values=np.asarray([pooled_loglik(lam,entries) for lam in grid])
    best=int(values.argmax());best_lambda=float(grid[best]);best_score=float(values[best])
    for i,lam in enumerate(grid):
        if (i>0 and values[i]<values[i-1]) or (i<len(grid)-1 and values[i]<values[i+1]):continue
        fit=minimize_scalar(lambda x:-pooled_loglik(x,entries),bounds=(max(0,lam-.01),min(1,lam+.01)),method='bounded',options={'xatol':1e-11})
        if -fit.fun>best_score:best_lambda=float(fit.x);best_score=float(-fit.fun)
    pools[assay]=dict(lambda_value=best_lambda,genes=len(entries),sum_profile_log_likelihood=best_score,
                      score_at_zero=float(values[0]),score_at_one=float(values[-1]),
                      method='Sum of gene-specific Gaussian profile log likelihoods; composite objective')
    profiles.extend(dict(assay=assay,lambda_value=float(lam),sum_profile_log_likelihood=float(v),relative_to_max=float(v-best_score)) for lam,v in zip(grid,values))
    assert best_score>=values.max()-1e-10
    for p in entries:
        f=fit_case(p,best_lambda);r=result_lookup[(assay,p['case']['gene'])]
        r.update(individual_lambda=r['lambda_ml'],individual_beta=r['beta'],individual_p=r['p'],individual_q=r['q'])
        beta=float(f['beta'][0,0]);se=float(f['se'][0,0]);df=int(f['df'][0]);crit=t.ppf(.975,df)
        r.update(lambda_ml=best_lambda,beta=beta,se=se,ci_low=float(beta-crit*se),ci_high=float(beta+crit*se),
                 t_stat=beta/se,df=df,p=float(f['p'][0]),q=None,intercept=float(f['intercept'][0]),
                 log_likelihood=float(f['log_likelihood'][0]),gls_r2=float(f['r_squared'][0]))

for assay in ['rna','atac']:
    fitted=[r for r in results if r['assay']==assay and r['p'] is not None]
    qs=false_discovery_control([r['p'] for r in fitted],method='bh')
    for r,q in zip(fitted,qs):r['q']=float(q)
    assert len({r['lambda_ml'] for r in fitted})==1
    payload['summary'][assay].update(
        nominal=[r['gene'] for r in fitted if r['p']<.05],
        fdr=[r['gene'] for r in fitted if r['q']<.05],
        prior_individual_lambda_fdr=prior['summary'][assay]['fdr'],
        shared_lambda=pools[assay]['lambda_value'])
for r in results:
    for k in ['individual_lambda','individual_beta','individual_p','individual_q']:r.setdefault(k,None)
payload['lambda_mode']='shared_by_assay'
payload['shared_lambda']=pools
payload['previous_analysis']=str(PREVIOUS/'analysis.json')
payload['previous_analysis_sha256']=hashlib.sha256((PREVIOUS/'analysis.json').read_bytes()).hexdigest()
payload['notes'][3]='One lambda is shared across all estimable genes within each assay. It is selected in [0,1] by maximizing the sum of full-model gene-specific profile log likelihoods, using a 0.01 grid and bounded refinement around every grid maximum; both endpoints are checked. Each gene retains its own intercept, group effect and residual variance.'
payload['notes'].insert(4,'The pooled objective is a composite likelihood: correlated genes and within-species rank dependencies are not modeled jointly. No confidence interval is assigned to the shared lambda. RNA and ATAC are separate assays with separate shared estimates.')
payload['notes'][6]='P-values and 95% intervals use the conventional GLS t approximation with n-2 degrees of freedom, conditional on the pooled lambda. They do not propagate uncertainty or cross-gene dependence in estimating lambda. Rank scores are bounded and tied, so these tests remain exploratory.'
payload['notes'].append('The same tree does not require equal Pagel lambda for different traits. Equal lambda across genes is imposed here at user request; the earlier gene-specific estimates are retained for comparison.')
assert payload['measurements']==prior['measurements']
assert payload['species']==prior['species'] and payload['panels']==prior['panels']
payload['references'].append('https://search.r-project.org/CRAN/refmans/caper/html/pgls.html')
(OUT/'analysis.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
validation=[]
for p in prepared:
    c=copy.deepcopy(p['case']);c['expected']=result_lookup[(c['assay'],c['gene'])];validation.append(c)
(OUT/'validation_cases.json').write_text(json.dumps(validation,indent=2,allow_nan=False)+'\n')
def write_tsv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=keys,delimiter='\t');w.writeheader()
        for r in rows:w.writerow({k:'NA' if r.get(k) is None else r[k] for k in keys})
write_tsv(OUT/'Shared_lambda_PGLS_all_genes.tsv',results)
write_tsv(OUT/'Shared_lambda_profiles.tsv',profiles)
write_tsv(OUT/'Shared_lambda_PGLS_measurements.tsv',payload['measurements'])
hits=[r for r in results if r['q'] is not None and r['q']<.05]
if hits:write_tsv(OUT/'Shared_lambda_PGLS_FDR05.tsv',hits)
print(json.dumps(dict(shared_lambda=pools,summary=payload['summary'],SIRT6=[r for r in results if r['gene']=='SIRT6']),indent=2),flush=True)

os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'outputs/sirt6_rearrangement_user_groups_2026-10-31/build/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
for ax,assay,title in zip(axes,['rna','atac'],['Liver RNA expression','TSS ±5 kb ATAC activity']):
    rows=[r for r in profiles if r['assay']==assay]
    ax.plot([r['lambda_value'] for r in rows],[r['relative_to_max'] for r in rows],c='#0072B2')
    ax.axvline(pools[assay]['lambda_value'],c='#D59B00',ls='--')
    ax.set_title(f'{title}\nShared λ = {pools[assay]["lambda_value"]:.4f}; {pools[assay]["genes"]} genes')
    ax.set_xlabel('Shared Pagel λ');ax.set_ylabel('Summed profile log likelihood\nrelative to maximum')
    ax.set_axisbelow(True);ax.grid(color='#E5E7EB',lw=.7)
fig.suptitle('One phylogenetic parameter across genes within each assay',fontsize=13)
for ext in ['png','svg','pdf']:fig.savefig(OUT/f'Shared_lambda_profiles.{ext}',dpi=220,facecolor='white')
plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(10,4.6),layout='constrained')
for ax,assay,title in zip(axes,['rna','atac'],['Liver RNA expression','TSS ±5 kb ATAC activity']):
    before=next(r for r in prior['results'] if r['assay']==assay and r['gene']=='SIRT6')
    after=result_lookup[(assay,'SIRT6')]
    for j,(r,color) in enumerate([(before,'#89949F'),(after,'#0072B2')]):
        ax.errorbar(r['beta'],j,xerr=r['ci_high']-r['beta'],fmt='o',c=color,capsize=4,ms=7)
    ax.set_yticks([0,1],['Gene-specific λ','Shared λ'])
    ax.set_ylim(1.6,-.6);ax.axvline(0,c='#333333',ls='--',lw=1)
    ax.set_title(f'{title}\nShared λ = {after["lambda_ml"]:.4f}\np = {after["p"]:.4g}; BH q = {after["q"]:.4g}')
    ax.set_xlabel('SIRT6 group effect and 95% interval\nPositive = higher score in rearranged species')
    ax.set_axisbelow(True);ax.xaxis.grid(color='#E5E7EB',lw=.7)
fig.suptitle('SIRT6 rank-score PGLS: shared versus gene-specific lambda',fontsize=13)
for ext in ['png','svg','pdf']:fig.savefig(OUT/f'SIRT6_shared_lambda_comparison.{ext}',dpi=220,facecolor='white')
plt.close(fig)
