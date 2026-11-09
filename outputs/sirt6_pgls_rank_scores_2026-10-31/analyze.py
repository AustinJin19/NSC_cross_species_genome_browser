"""PGLS of within-species rank score on rearrangement group, with per-gene ML lambda."""
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import t
from Bio import Phylo

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SOURCE = ROOT/'outputs/sirt6_rearrangement_user_groups_2026-10-31/analysis.json'
TREE_MATRIX = ROOT/'cross_species_visual_tool/web_app/data/phylogenetic_covariance.json'
TREE_FILE = ROOT/'cross_species_visual_tool/web_app/data/timetree_pgls.nwk'
sys.path.insert(0, str(ROOT/'outputs/pgls_genomewide_2026-10-27'))
from pgls_batch import fit_batch

(OUT/'build').mkdir(exist_ok=True)
source = json.loads(SOURCE.read_text())
tree = json.loads(TREE_MATRIX.read_text())
species = {s['id']:s for s in source['species']}
assert source['group_counts'] == {'Rearranged':10, 'Unrearranged':20, 'Excluded':1}
assert all(s['latin'] in tree['species'] for s in species.values())
tree_index = {name:i for i,name in enumerate(tree['species'])}
full_C = np.asarray(tree['covariance'])
observations = defaultdict(list)
measurements = []
for r in source['measurements']:
    row = dict(r, latin=species[r['species_id']]['latin'], rank_score=None,
               rearranged=int(r['group']=='Rearranged') if r['group']!='Excluded' else None)
    if r['rank'] is not None and r['reference_genes'] >= 2:
        row['rank_score'] = 1-(r['rank']-1)/(r['reference_genes']-1)
        assert 0 <= row['rank_score'] <= 1
    if r['used_in_test']:
        assert row['rank_score'] is not None
        observations[(r['assay'],r['gene'])].append(row)
    measurements.append(row)

results = []
validation_cases = []
for prior in source['results']:
    key=(prior['assay'],prior['gene'])
    rows=sorted(observations[key],key=lambda x:x['latin'])
    result={k:prior[k] for k in ['assay','gene','reference_genes','n_rearranged','n_unrearranged']}
    result.update(status=prior['status'],n=len(rows),lambda_ml=None,beta=None,se=None,ci_low=None,ci_high=None,
                  t_stat=None,df=None,p=None,q=None,intercept=None,log_likelihood=None,gls_r2=None,
                  score_median_rearranged=None,score_median_unrearranged=None,unique_scores=None,
                  ols_beta=None,ols_p=None,bm_beta=None,bm_p=None,
                  wilcoxon_p=prior['p'],wilcoxon_q=prior['q'])
    if prior['status']!='Tested':
        results.append(result)
        continue
    y=np.asarray([r['rank_score'] for r in rows])
    x=np.asarray([r['rearranged'] for r in rows])
    result.update(unique_scores=len(np.unique(y)),score_median_rearranged=float(np.median(y[x==1])),
                  score_median_unrearranged=float(np.median(y[x==0])))
    if np.ptp(y)<1e-14:
        result['status']='Not estimable: constant rank score'
        results.append(result)
        continue
    indexes=[tree_index[r['latin']] for r in rows]
    C=full_C[np.ix_(indexes,indexes)]
    X=np.column_stack([np.ones(len(y)),x])
    fitted={mode:fit_batch(y[None,:],X[None,:,:],C,mode) for mode in ['ML',0,1]}
    f=fitted['ML']
    if not f['valid'][0]:
        result['status']='Not estimable: singular or zero residual variance'
        results.append(result)
        continue
    beta=float(f['beta'][0,0]);se=float(f['se'][0,0]);df=int(f['df'][0]);crit=t.ppf(.975,df)
    result.update(status='Tested',lambda_ml=float(f['lambda_value'][0]),beta=beta,se=se,
                  ci_low=float(beta-crit*se),ci_high=float(beta+crit*se),t_stat=beta/se,df=df,
                  p=float(f['p'][0]),intercept=float(f['intercept'][0]),
                  log_likelihood=float(f['log_likelihood'][0]),gls_r2=float(f['r_squared'][0]),
                  ols_beta=float(fitted[0]['beta'][0,0]),ols_p=float(fitted[0]['p'][0]),
                  bm_beta=float(fitted[1]['beta'][0,0]),bm_p=float(fitted[1]['p'][0]))
    # Check the lambda=0 fit against ordinary least squares independently.
    beta_ols=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@beta_ols
    se_ols=np.sqrt(np.diag(np.linalg.inv(X.T@X))*np.dot(resid,resid)/(len(y)-2))
    assert abs(beta_ols[1]-result['ols_beta'])<1e-10
    assert abs(2*t.sf(abs(beta_ols[1]/se_ols[1]),len(y)-2)-result['ols_p'])<1e-9
    validation_cases.append(dict(assay=prior['assay'],gene=prior['gene'],points=[dict(species=r['latin'].replace(' ','_'),x=r['rearranged'],y=r['rank_score']) for r in rows],
                                 covariance=C.tolist(),expected=result.copy()))
    results.append(result)

def bh(rows):
    ordered=sorted([r for r in rows if r['p'] is not None],key=lambda r:r['p'])
    q=1
    for i in range(len(ordered)-1,-1,-1):
        q=min(q,ordered[i]['p']*len(ordered)/(i+1));ordered[i]['q']=q
for assay in ['rna','atac']:
    bh([r for r in results if r['assay']==assay])
summary={a:dict(tested=sum(r['status']=='Tested' for r in results if r['assay']==a),
                constant=sum(r['status']=='Not estimable: constant rank score' for r in results if r['assay']==a),
                nominal=[r['gene'] for r in results if r['assay']==a and r['p'] is not None and r['p']<.05],
                fdr=[r['gene'] for r in results if r['assay']==a and r['q'] is not None and r['q']<.05],
                prior_wilcoxon_fdr=source['summary'][a]['fdr05']) for a in ['rna','atac']}

# Inspect independent phylogenetic replication of group labels.
phy=Phylo.read(TREE_FILE,'newick')
included_names={s['latin'].replace(' ','_') for s in species.values() if s['group']!='Excluded'}
rearranged_names={s['latin'].replace(' ','_') for s in species.values() if s['group']=='Rearranged'}
mrca=phy.common_ancestor(sorted(rearranged_names))
descendants={x.name for x in mrca.get_terminals()} & included_names
clade=dict(rearranged_monophyletic_in_included_tree=descendants==rearranged_names,
           included_descendants=sorted(descendants),included_unrearranged_within_clade=sorted(descendants-rearranged_names))

notes=[
 'Outcome: rank score = 1 - (shared-panel rank - 1)/(shared-panel gene count - 1). A score of 1 is highest and 0 is lowest. Ties retain average ranks.',
 'The shared reference panel, missing values and species-median measurements are unchanged from the user-group Wilcoxon analysis. Panels are target-specific; each test uses the same reference genes across its included species.',
 'Model: rank score = intercept + beta * rearranged + error; rearranged=1, unrearranged=0. Positive beta means higher rank score in rearranged species. No MLS, body mass, length or batch predictor is included.',
 'Pagel lambda is estimated independently for each gene and assay by Gaussian maximum likelihood in [0,1], checking endpoints and local maxima. The covariance is the corresponding submatrix of the existing TimeTree 6 matrix, preserving its full-tree root.',
 'Covariance diagonal = root-to-tip branch length; off diagonal = lambda times shared root-to-MRCA length. Branch lengths are million years. Lambda describes residual dependence conditional on group, not dependence of the unadjusted outcome.',
 'P-values and 95% intervals use the conventional GLS t approximation with n-2 degrees of freedom and residual variance SSE/(n-2). They do not integrate uncertainty in estimated lambda. Rank scores are bounded and tied; Gaussian inference is exploratory, especially with few unique scores.',
 'BH FDR is calculated separately for valid RNA and ATAC PGLS fits. Constant-score genes have no estimable Gaussian association test and are excluded from the FDR family. The prior Wilcoxon family included constant ranks with p=1, so family sizes can differ.',
 'Species assignments: 10 rearranged, 20 unrearranged. Muskrat, human, macaque and bowhead whale are excluded. The seven user-resolved marker assignments are preserved and are not treated as independent annotation verification.',
 'All ten rearranged species form one clade among included taxa. This comparison may represent one inherited rearrangement, not ten independent events; PGLS cannot distinguish rearrangement causation from other changes on that ancestral lineage.' if clade['rearranged_monophyletic_in_included_tree'] else 'Shared ancestry may reduce the number of independent rearrangement events; PGLS does not prove causation.',
 'SIRT6 RNA ranks use only three shared genes (SIRT6, CHAF1A, HDGFL2) and have coarse tied scores. SIRT6 ATAC ranks use 30 shared genes. Compare assays with this coverage difference in mind.',
 'RNA uses source liver log.tpm; ATAC is integrated TSS ±5 kb CPM·bp before ranking. Exact original labels and highest-median loci are retained. The source available-window rule is at least 50%; missing values are never imputed as zero.',
 'These are saved external statistical results. Rerun analyze.py and the workbook builder after changing data, assignments or scoring. Lambda=0 (OLS) and lambda=1 (Brownian) fits are sensitivity checks, not additional primary discoveries.',
]
payload=dict(results=results,summary=summary,measurements=measurements,species=source['species'],
             panels=source['panels'],notes=notes,group_counts=source['group_counts'],clade=clade,
             source=str(SOURCE),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             measurement_source=source['source'],measurement_sha256=source['source_sha256'],
             tree_file=str(TREE_FILE),tree_matrix=str(TREE_MATRIX),tree_sha256=hashlib.sha256(TREE_MATRIX.read_bytes()).hexdigest(),
             analysis_date='2026-10-31',references=[tree['source_url'],'https://search.r-project.org/CRAN/refmans/phylolm/html/phylolm.html','https://search.r-project.org/CRAN/refmans/ape/html/corPagel.html'])
(OUT/'analysis.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
(OUT/'validation_cases.json').write_text(json.dumps(validation_cases,indent=2,allow_nan=False)+'\n')
def write_tsv(path,rows):
    with path.open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');w.writeheader()
        for r in rows:w.writerow({k:'NA' if v is None else v for k,v in r.items()})
write_tsv(OUT/'PGLS_rank_score_all_genes.tsv',results)
write_tsv(OUT/'PGLS_rank_score_measurements.tsv',measurements)
hits=[r for r in results if r['q'] is not None and r['q']<.05]
if hits:write_tsv(OUT/'PGLS_rank_score_FDR05.tsv',hits)
print(json.dumps(dict(summary=summary,SIRT6=[r for r in results if r['gene']=='SIRT6'],clade=clade),indent=2),flush=True)

os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'outputs/sirt6_rearrangement_user_groups_2026-10-31/build/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
colors={'Rearranged':'#0072B2','Unrearranged':'#D59B00'}
fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
for ax,assay,title in zip(axes,['rna','atac'],['Liver RNA expression','TSS ±5 kb ATAC activity']):
    f=next(r for r in results if r['assay']==assay and r['gene']=='SIRT6')
    obs=observations[(assay,'SIRT6')]
    for j,group in enumerate(['Unrearranged','Rearranged']):
        values=[r['rank_score'] for r in obs if r['group']==group]
        ax.boxplot(values,positions=[j],widths=.4,patch_artist=True,showfliers=False,
            boxprops={'facecolor':colors[group]+'25','edgecolor':colors[group]},medianprops={'color':'#222222'},
            whiskerprops={'color':colors[group]},capprops={'color':colors[group]})
        offsets=np.random.default_rng(12+j).uniform(-.12,.12,len(values))
        ax.scatter(j+offsets,values,c=colors[group],s=35,edgecolor='white',linewidth=.5,zorder=3)
    ax.plot([0,1],[f['intercept'],f['intercept']+f['beta']],color='#262626',ls='--',lw=1.5,label='PGLS fitted group score')
    ax.set_xticks([0,1],[f'Unrearranged\n(n={f["n_unrearranged"]})',f'Rearranged\n(n={f["n_rearranged"]})'])
    ax.set_ylim(-.06,1.06);ax.set_ylabel(f'SIRT6 rank score ({f["reference_genes"]} shared genes)\n1 = highest, 0 = lowest')
    ax.set_title(f'{title}\nλ = {f["lambda_ml"]:.3f}; β = {f["beta"]:.3f}\np = {f["p"]:.4g}; BH q = {f["q"]:.4g}',pad=12)
    ax.set_axisbelow(True);ax.yaxis.grid(True,color='#E5E7EB',lw=.7);ax.legend(loc='lower left',fontsize=8)
fig.suptitle('SIRT6 rank scores with phylogenetic regression',fontsize=14)
for ext in ['png','svg','pdf']:fig.savefig(OUT/f'SIRT6_PGLS_rank_scores.{ext}',dpi=220,facecolor='white')
plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,10),layout='constrained')
for ax,assay,title in zip(axes,['rna','atac'],['Liver RNA expression','TSS ±5 kb ATAC activity']):
    rows=sorted([r for r in results if r['assay']==assay and r['p'] is not None],key=lambda r:r['beta'])
    y=np.arange(len(rows));bs=np.asarray([r['beta'] for r in rows]);err=np.asarray([1.0*(r['ci_high']-r['beta']) for r in rows])
    ax.errorbar(bs,y,xerr=err,fmt='none',ecolor='#A5ACB5',lw=.8,zorder=1)
    ax.scatter(bs,y,c=['#0072B2' if r['q']<.05 else '#9EA7B0' for r in rows],s=24,zorder=2)
    ax.set_yticks(y,[r['gene']+(' **' if r['q']<.05 else ' *' if r['p']<.05 else '') for r in rows],fontsize=8)
    for tick,r in zip(ax.get_yticklabels(),rows):
        if r['gene']=='SIRT6':tick.set_fontweight('bold');tick.set_color('#A51624')
    ax.axvline(0,c='#222222',lw=.8);ax.set_title(title)
    ax.set_xlabel('PGLS group coefficient and 95% interval\nPositive = higher rank score in rearranged species')
    ax.set_axisbelow(True);ax.xaxis.grid(True,color='#E5E7EB',lw=.7)
fig.suptitle('Phylogenetically adjusted rank-score differences\n* nominal p < 0.05; ** BH q < 0.05',fontsize=13)
for ext in ['png','svg','pdf']:fig.savefig(OUT/f'Locus_PGLS_rank_score_effects.{ext}',dpi=220,facecolor='white')
plt.close(fig)
