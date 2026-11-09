"""Within-group concordance of SIRT6-locus rankings, with missingness explicit."""
import csv
import hashlib
import itertools
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist
from scipy.stats import chi2, false_discovery_control, friedmanchisquare, rankdata, spearmanr

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
(OUT/'build').mkdir(exist_ok=True)
RNA_FILE=ROOT/'outputs/sirt6_rna_tissues_2026-10-31/analysis.json'
ATAC_FILE=ROOT/'outputs/sirt6_shared_lambda_rank_scores_2026-10-31/analysis.json'
rna=json.loads(RNA_FILE.read_text());atac=json.loads(ATAC_FILE.read_text())
genes=rna['genes'];species=[s for s in rna['species'] if s['group']!='Excluded']
species_by_id={s['id']:s for s in species}
GROUPS=['Rearranged','Unrearranged']
datasets=[]
for tissue in rna['tissues']:
    records=[r for r in rna['measurements'] if r['tissue']==tissue]
    values={(r['species_id'],r['gene']):r['median_log_tpm'] for r in records
        if r['species_id'] in species_by_id and r['median_log_tpm'] is not None}
    datasets.append(dict(dataset='RNA '+tissue,assay='RNA',tissue=tissue,values=values,expected_panel=rna['panels'][tissue+'|SIRT6']))
values={(r['species_id'],r['gene']):r['value'] for r in atac['measurements']
    if r['assay']=='atac' and r['species_id'] in species_by_id and r['value'] is not None}
datasets.append(dict(dataset='ATAC TSS ±5 kb',assay='ATAC',tissue='Liver',values=values,expected_panel=atac['panels']['atac|SIRT6']))

def concordance(R):
    m,k=R.shape
    ties=sum(sum(c**3-c for c in Counter(row).values()) for row in R)
    correction=1-ties/(m*(k**3-k))
    if correction<=0:return None
    S=float(np.sum((R.sum(axis=0)-m*(k+1)/2)**2))
    W=12*S/(m*m*(k**3-k)-m*ties)
    Q=m*(k-1)*W
    scipy=friedmanchisquare(*R.T)
    assert abs(Q-scipy.statistic)<1e-10
    assert 0<=W<=1+1e-12
    correlations=[float(spearmanr(a,b).statistic) for a,b in itertools.combinations(R,2)]
    return dict(kendall_w=min(1.,W),friedman_chi2=Q,df=k-1,p_chi2=float(chi2.sf(Q,k-1)),
        tie_correction=correction,mean_pairwise_rho=float(np.mean(correlations)),
        min_pairwise_rho=float(min(correlations)),max_pairwise_rho=float(max(correlations)))

def exact_probability(R):
    """Exact upper-tail rank-label randomization by DP on rank sums (K <= 3).

    Fixing the first row removes redundant global column permutations.
    Duplicate tied permutations retain their multiplicity.
    """
    m,k=R.shape
    doubled=np.rint(2*R).astype(int)
    observed=int(np.sum(doubled.sum(axis=0)**2))
    states={tuple(doubled[0,:-1]):1}
    for row in doubled[1:]:
        possible=Counter(itertools.permutations(row.tolist()))
        updated=defaultdict(int)
        for state,count in states.items():
            for perm,multiplicity in possible.items():
                updated[tuple(a+b for a,b in zip(state,perm[:-1]))]+=count*multiplicity
        states=updated
    total=math.factorial(k)**(m-1)
    assert sum(states.values())==total
    total_rank=m*k*(k+1)
    exceed=sum(count for state,count in states.items()
        if sum(x*x for x in state)+(total_rank-sum(state))**2>=observed)
    return dict(p_permutation=exceed/total,permutation_method='Exact rank-label permutations',
        permutations=total,exceedances=exceed,mc_upper95=None,seed=None)

def monte_carlo_probability(R,seed,B=200000):
    rng=np.random.default_rng(seed)
    observed=np.sum(R.sum(axis=0)**2)
    count=0
    for start in range(0,B,1000):
        size=min(1000,B-start)
        # Each species independently permutes its gene labels; ties are retained.
        permutations=rng.permuted(np.broadcast_to(R,(size,*R.shape)).copy(),axis=2)
        statistic=np.sum(permutations.sum(axis=1)**2,axis=1)
        count+=int(np.sum(statistic>=observed-1e-9))
    return dict(p_permutation=(count+1)/(B+1),permutation_method='Monte Carlo rank-label permutations',
        permutations=B,exceedances=count,mc_upper95=float(beta_dist.ppf(.95,count+1,B-count)) if count<B else 1.,seed=seed)

# Validate the exact DP against direct enumeration on small tied examples.
for example in [np.array([[1,2,3],[1,3,2],[2,1,3]]),np.array([[1.5,1.5,3],[1,2.5,2.5],[3,1.5,1.5]])]:
    target=np.sum(example.sum(axis=0)**2)
    possibilities=list(itertools.product(*[list(itertools.permutations(r)) for r in example[1:]]))
    direct=sum(np.sum(np.vstack([example[0],*p]).sum(axis=0)**2)>=target-1e-10 for p in possibilities)/len(possibilities)
    assert abs(direct-exact_probability(example)['p_permutation'])<1e-12

summary=[];profiles=[];pairs=[];coverage=[];cases=[];panels={};consensus=[]
for di,ds in enumerate(datasets):
    key=ds['dataset'];values=ds['values']
    available=[s for s in species if any((s['id'],g) in values for g in genes)]
    eligible=[s for s in species if (s['id'],'SIRT6') in values]
    common=[g for g in genes if eligible and all((s['id'],g) in values for s in eligible)]
    assert common==ds['expected_panel']
    panels[key]=common
    for s in species:
        coverage.append(dict(dataset=key,species=s['common'],species_id=s['id'],group=s['group'],
            measured_genes=sum((s['id'],g) in values for g in genes),sirt6_measured=s in eligible,
            primary_profile_included=s in eligible,reason='Included in SIRT6 comparison panel' if s in eligible else 'SIRT6 unavailable'))
    for group in GROUPS:
        ss=[s for s in eligible if s['group']==group]
        all_group=[s for s in available if s['group']==group]
        group_all_common=[g for g in genes if all_group and all((s['id'],g) in values for s in all_group)]
        r=dict(dataset=key,assay=ds['assay'],tissue=ds['tissue'],group=group,n_species=len(ss),n_genes=len(common),
            all_measured_species=len(all_group),genes_shared_by_all_measured=len(group_all_common),status='Tested',
            kendall_w=None,friedman_chi2=None,df=None,p_chi2=None,q_chi2=None,p_permutation=None,q_permutation=None,
            permutation_method=None,permutations=None,exceedances=None,mc_upper95=None,seed=None,
            tie_correction=None,mean_pairwise_rho=None,min_pairwise_rho=None,max_pairwise_rho=None)
        if ss and common:
            Y=np.array([[values[s['id'],g] for g in common] for s in ss])
            R=rankdata(-Y,axis=1,method='average')
            for i,s in enumerate(ss):
                for j,g in enumerate(common):
                    profiles.append(dict(dataset=key,group=group,species=s['common'],species_id=s['id'],gene=g,
                        source_median=float(Y[i,j]),rank=float(R[i,j]),rank_score=float(1-(R[i,j]-1)/(len(common)-1)) if len(common)>1 else None))
            for j,g in enumerate(common):
                consensus.append(dict(dataset=key,group=group,gene=g,mean_rank=float(np.mean(R[:,j])),
                    median_rank=float(np.median(R[:,j])),min_rank=float(np.min(R[:,j])),max_rank=float(np.max(R[:,j]))))
        if len(ss)<3:
            r['status']='Not tested: fewer than 3 species'
        elif len(common)<3:
            r['status']='Not tested: fewer than 3 shared genes'
        else:
            c=concordance(R)
            if c is None:r['status']='Not estimable: all genes tied'
            else:
                r.update(c)
                r.update(exact_probability(R) if len(common)<=3 else monte_carlo_probability(R,20261031+di*2+GROUPS.index(group)))
                cases.append(dict(dataset=key,group=group,ranks=R.tolist(),genes=common,species=[s['latin'] for s in ss],expected=r))
        summary.append(r)
    for a,b in itertools.combinations(available,2):
        shared=[g for g in genes if (a['id'],g) in values and (b['id'],g) in values]
        rho=None;status='Tested descriptively'
        if len(shared)<3:status='Fewer than 3 shared genes'
        else:
            x=np.array([values[a['id'],g] for g in shared]);y=np.array([values[b['id'],g] for g in shared])
            if np.ptp(x)==0 or np.ptp(y)==0:status='Constant values'
            else:rho=float(spearmanr(x,y).statistic)
        pairs.append(dict(dataset=key,species_a=a['common'],species_b=b['common'],species_id_a=a['id'],species_id_b=b['id'],
            group_a=a['group'],group_b=b['group'],comparison=a['group'] if a['group']==b['group'] else 'Between groups',
            shared_genes=len(shared),spearman_rho=rho,status=status,gene_names='; '.join(shared)))
    print(key,[(r['group'],r['n_species'],r['n_genes'],r['kendall_w'],r['p_permutation']) for r in summary if r['dataset']==key],flush=True)

valid=[r for r in summary if r['status']=='Tested']
for p,q in [('p_permutation','q_permutation'),('p_chi2','q_chi2')]:
    for r,fdr in zip(valid,false_discovery_control([r[p] for r in valid],method='bh')):r[q]=float(fdr)
pair_summary=[]
for ds in datasets:
    for comparison in GROUPS+['Between groups']:
        rr=[r for r in pairs if r['dataset']==ds['dataset'] and r['comparison']==comparison and r['spearman_rho'] is not None]
        pair_summary.append(dict(dataset=ds['dataset'],comparison=comparison,pairs=len(rr),
            median_rho=float(np.median([r['spearman_rho'] for r in rr])) if rr else None,
            median_shared_genes=float(np.median([r['shared_genes'] for r in rr])) if rr else None,
            min_shared_genes=min(r['shared_genes'] for r in rr) if rr else None,max_shared_genes=max(r['shared_genes'] for r in rr) if rr else None))

notes=[
 'Question: do species within each rearrangement group share a similar ordering of locus genes? This is a within-group concordance analysis, not a gene-wise difference, variance or MLS test.',
 'Primary panels preserve the prior SIRT6 comparison: only species with a SIRT6 measurement are eligible. Within each dataset, both groups rank exactly the same genes measured in all those species. RNA uses 3 genes; ATAC TSS ±5 kb uses 30. Species inclusion and genes are listed explicitly.',
 'Across all RNA species with any locus measurement, there is no single gene measured in every unrearranged species. A complete-block Friedman analysis of that full set is not possible without discarding or imputing data. The SIRT6 panel is a restricted, previously defined subset; no missing values are imputed.',
 'Rank 1 denotes highest species-median RNA log.tpm or ATAC integrated CPM·bp. Ties receive average ranks. One species contributes one profile, irrespective of sample count.',
 'Kendall W measures concordance of gene ordering across species, with 1 denoting identical order. Species are blocks/raters and genes are treatments/objects. Ties are corrected. Friedman chi-square Q = M*(K-1)*W, with df=K-1.',
 'A small Friedman p-value indicates a systematic gene hierarchy across species; it does not prove all species have identical rankings. Agreement strength is described by W. Significant concordance in both groups does not establish that the groups share the same hierarchy.',
 'Because RNA panels are small, the primary p-value is the exact upper-tail probability under independent gene-label permutations within species. A dynamic program sums all rank-sum configurations, retaining tie multiplicities; one row is fixed by global label symmetry.',
 'For the 30-gene ATAC panel, use 200,000 independent rank-label randomizations with recorded seeds and p=(exceedances+1)/(B+1). No simulated p-value is reported as zero. The Monte Carlo upper95 field is an upper 95% binomial bound on the tail probability.',
 'BH correction pools all valid dataset-by-group concordance tests. The asymptotic Friedman chi-square p and its separate BH q are supplied for reference; permutation p is primary.',
 'Pairwise Spearman correlations use all measured locus genes shared by each species pair (minimum 3), with counts and gene names recorded. Gene sets vary across pairs. Pairwise correlations are descriptive, not independent replicates, and have no pooled significance test.',
 'These concordance tests are not corrected for phylogeny and do not estimate lambda. Related species can share rank patterns through ancestry. All rearranged species form one clade; p-values are exploratory under the independent-species null.',
 'Human, macaque, bowhead whale and muskrat remain excluded. Tissue correction BMR_Turk_Lu to lung, source priority, species medians, symbol matching and user group assignments are unchanged from the prior analyses.',
 'No conclusion about all 63 genes should be drawn from the three-gene RNA panel (SIRT6, CHAF1A, HDGFL2). Larger pairwise overlaps are provided separately and do not repair missing complete profiles.',
]
payload=dict(summary=summary,pair_summary=pair_summary,profiles=profiles,pairs=pairs,coverage=coverage,consensus=consensus,panels=panels,
    species=species,datasets=[d['dataset'] for d in datasets],notes=notes,
    sources=[dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [RNA_FILE,ATAC_FILE]],
    references=['https://stat.ethz.ch/R-manual/R-devel/library/stats/html/friedman.test.html',
        'https://search.r-project.org/CRAN/refmans/DescTools/html/KendallW.html',
        'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html'])
(OUT/'analysis.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
(OUT/'validation_cases.json').write_text(json.dumps(cases,indent=2,allow_nan=False)+'\n')
def tsv(name,rr):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader()
        w.writerows({k:'NA' if v is None else v for k,v in r.items()} for r in rr)
for name,rr in [('Concordance_tests',summary),('Pairwise_summary',pair_summary),('Profile_ranks',profiles),('Pairwise_correlations',pairs),('Species_coverage',coverage),('Consensus_gene_ranks',consensus)]:tsv(name+'.tsv',rr)

os.environ.setdefault('MPLCONFIGDIR',str(OUT/'build/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
shown=[d['dataset'] for d in datasets if any(r['dataset']==d['dataset'] and r['status']=='Tested' for r in summary)]
fig,ax=plt.subplots(figsize=(9,5.4),layout='constrained')
colors={'Rearranged':'#0072B2','Unrearranged':'#D59B00'}
for gi,group in enumerate(GROUPS):
    rr=[next(r for r in summary if r['dataset']==ds and r['group']==group) for ds in shown]
    positions=np.arange(len(shown))+(gi-.5)*.23
    ax.scatter([r['kendall_w'] for r in rr],positions,c=colors[group],label=group,s=60,zorder=3)
    for p,r in zip(positions,rr):
        ax.text(r['kendall_w']-.025,p,f"{r['kendall_w']:.3f}",ha='right',va='center',fontsize=8,color=colors[group])
ax.set_yticks(range(len(shown)),[ds+(' (30 genes)' if ds.startswith('ATAC') else ' (3 genes)') for ds in shown])
ax.invert_yaxis();ax.set_xlim(0,1.05);ax.set_xlabel("Kendall’s W: higher = stronger agreement in gene order")
ax.set_title('Within-group concordance of locus gene rankings\nSIRT6-comparable species; no phylogenetic correction',fontsize=12)
ax.set_axisbelow(True);ax.xaxis.grid(color='#E5E7EB');ax.legend(loc='lower left',frameon=False)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'Within_group_rank_concordance.{ext}',dpi=220,facecolor='white',bbox_inches='tight')
plt.close(fig)
