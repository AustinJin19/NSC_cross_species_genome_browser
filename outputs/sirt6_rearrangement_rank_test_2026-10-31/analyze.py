"""Gene-order association with within-species ranks; no browser data are modified."""
import csv
import argparse
import hashlib
import itertools
import json
import math
import os
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, rankdata

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--out', type=Path, default=SCRIPT_DIR)
parser.add_argument('--assignments', type=Path)
args = parser.parse_args()
OUT = args.out.resolve()
(OUT/'build').mkdir(parents=True, exist_ok=True)
overrides = json.loads(args.assignments.read_text()) if args.assignments else {}
SOURCE = ROOT/'outputs/sirt6_neighborhood_rank_matrices_2026-10-29/matrix_data.json'
sys.path.insert(0, str(ROOT/'cross_species_visual_tool/web_app'))
import server

EXCLUDED = {'Human', 'Macaque', 'Balaena_mysticetus'}
MARKERS = ['SIRT6', 'CREB3L3', 'TLE6']
GROUPS = ['Rearranged', 'Unrearranged']
server.load_catalog()
data = json.loads(SOURCE.read_text())
species = [dict(s, index=i) for i, s in enumerate(data['metadata']['species']) if s['id'] not in EXCLUDED]
genes = data['metadata']['genes']
assert len(species) == 31 and len(genes) == 63

classification = []
for sp in species:
    record = {k:sp[k] for k in ['id','common','latin']}
    record.update(group='Unclassified', reason='', evidence='Local annotation gene order; provisional')
    loci = {g:server.GENES[sp['id']].get(g, []) for g in MARKERS}
    for g in MARKERS:
        record[g+'_loci'] = loci[g]
    if any(len(v) != 1 for v in loci.values()):
        record['reason'] = '; '.join(f'{g}: {len(loci[g])} annotated loci' for g in MARKERS if len(loci[g]) != 1)
    else:
        a,b,c = [loci[g][0] for g in MARKERS]
        if len({g['chrom'] for g in [a,b,c]}) != 1:
            record['reason'] = 'Markers on different scaffolds/chromosomes; order unresolved'
        else:
            orientation = 1 if a['strand']=='+' else -1
            bpos,cpos = [(g['tss']-a['tss'])*orientation for g in [b,c]]
            record.update(CREB3L3_relative_bp=bpos, TLE6_relative_bp=cpos)
            if min(bpos,cpos)<=0 or bpos==cpos:
                record['reason'] = 'Different marker configuration; does not match either comparison group'
            else:
                record['group'] = 'Rearranged' if cpos<bpos else 'Unrearranged'
                record['reason'] = 'SIRT6 → TLE6 → CREB3L3' if cpos<bpos else 'SIRT6 → CREB3L3 → TLE6'
    record['annotation_group'] = record['group']
    record['annotation_reason'] = record['reason']
    if sp['id'] in overrides:
        override = overrides[sp['id']]
        assert override['group'] in GROUPS + ['Excluded']
        record.update(override)
    classification.append(record)
by_species = {r['id']:r for r in classification}
assert set(overrides) <= set(by_species), 'Unknown species in assignment overrides'


@lru_cache(None)
def sum_distribution(rank_twice, n):
    """Exact label-permutation distribution; duplicated ranks retain their multiplicity."""
    dp = [defaultdict(int) for _ in range(n+1)]
    dp[0][0] = 1
    for i, score in enumerate(rank_twice):
        for size in range(min(n, i+1), 0, -1):
            for previous, count in list(dp[size-1].items()):
                dp[size][previous+score] += count
    assert sum(dp[n].values()) == math.comb(len(rank_twice), n)
    return dict(dp[n])


def exact_wilcoxon(x,y):
    """Two-sided Wilcoxon rank-sum/Mann–Whitney test, exact even with ties."""
    x,y = np.asarray(x,float),np.asarray(y,float)
    ranks = rankdata(np.concatenate([x,y]), method='average')
    twice = tuple(int(round(v*2)) for v in ranks)
    observed = sum(twice[:len(x)])
    distribution = sum_distribution(tuple(sorted(twice)),len(x))
    total = sum(distribution.values())
    lower = sum(n for score,n in distribution.items() if score<=observed)/total
    upper = sum(n for score,n in distribution.items() if score>=observed)/total
    u = observed/2-len(x)*(len(x)+1)/2
    return {'u':float(u),'p':min(1.,2*min(lower,upper)),
            'delta_higher':float(1-2*u/(len(x)*len(y))), 'permutations':total}


def bh(values):
    result = [None]*len(values)
    valid = sorted((p,i) for i,p in enumerate(values) if p is not None)
    previous=1.
    for rank,(p,i) in reversed(list(enumerate(valid,1))):
        previous=min(previous,p*len(valid)/rank)
        result[i]=previous
    return result


# Independent checks: enumeration with ties and SciPy's exact no-tie distribution.
for x,y in [([1,2,2],[2,4,5]),([2,2],[2,2,3]),([1,3,4],[2,5,6]),([1,1],[1,1])]:
    r=rankdata(x+y); observed=sum(r[:len(x)])
    sums=[sum(r[list(c)]) for c in itertools.combinations(range(len(r)),len(x))]
    p=min(1,2*min(sum(s<=observed for s in sums),sum(s>=observed for s in sums))/len(sums))
    assert abs(exact_wilcoxon(x,y)['p']-p)<1e-12
for n,m in [(3,3),(9,14),(6,8)]:
    values=np.random.default_rng(n+m).permutation(n+m).tolist()
    assert abs(exact_wilcoxon(values[:n],values[n:])['p']-mannwhitneyu(values[:n],values[n:],method='exact').pvalue)<1e-10

results=[]
measurements=[]
panels={}
for assay in ['rna','atac']:
    rows={r['gene']:r for r in data[assay]}
    values={sp['id']:{g:rows[g]['values'][sp['index']] for g in genes} for sp in species}
    original={}
    for sp in species:
        sid=sp['id']; valid={g:v for g,v in values[sid].items() if v is not None}
        ranks=rankdata([-v for v in valid.values()],method='average')
        original[sid]={g:float(rank) for g,rank in zip(valid,ranks)}
    for gene in genes:
        included=[sp for sp in species if by_species[sp['id']]['group'] in GROUPS and values[sp['id']][gene] is not None]
        common=[g for g in genes if included and all(values[sp['id']][g] is not None for sp in included)]
        panels[assay+'|'+gene]=common
        common_ranks={}
        for sp in included:
            sid=sp['id']
            common_ranks[sid]=float(rankdata([-values[sid][g] for g in common],method='average')[common.index(gene)])
        grouped={group:[sp['id'] for sp in included if by_species[sp['id']]['group']==group] for group in GROUPS}
        rec={'assay':assay,'gene':gene,'reference_genes':len(common),
             'n_rearranged':len(grouped['Rearranged']),'n_unrearranged':len(grouped['Unrearranged']),
             'p':None,'q':None,'status':'Tested','u':None,'delta_higher':None,'permutations':None,
             'median_rearranged':None,'median_unrearranged':None,
             'original_rank_p':None,'percentile_p':None,'raw_value_p':None}
        if min(len(v) for v in grouped.values())<3:
            rec['status']='Not tested: fewer than 3 measured species in a group'
        elif len(common)<2:
            rec['status']='Not tested: fewer than 2 shared genes for ranking'
        else:
            x=[common_ranks[sid] for sid in grouped['Rearranged']]
            y=[common_ranks[sid] for sid in grouped['Unrearranged']]
            rec.update(exact_wilcoxon(x,y))
            rec.update(median_rearranged=float(np.median(x)),median_unrearranged=float(np.median(y)))
            for key in ['original_rank','percentile','raw_value']:
                def score(sid):
                    if key=='raw_value':return values[sid][gene]
                    rank=original[sid][gene]
                    return rank if key=='original_rank' else (rank-1)/(len(original[sid])-1)
                rec[key+'_p']=exact_wilcoxon([score(sid) for sid in grouped['Rearranged']], [score(sid) for sid in grouped['Unrearranged']])['p']
        results.append(rec)
        for sp in species:
            sid=sp['id']; value=values[sid][gene]; oldrank=original[sid].get(gene)
            measurements.append(dict(assay=assay,gene=gene,species_id=sid,species=sp['common'],
                group=by_species[sid]['group'],value=value,original_rank=oldrank,
                original_genes_measured=len(original[sid]),rank=common_ranks.get(sid),reference_genes=len(common),
                used_in_test=rec['status']=='Tested' and sid in common_ranks))
    current=[r for r in results if r['assay']==assay]
    for key in ['p','original_rank_p','percentile_p','raw_value_p']:
        for r,q in zip(current,bh([r[key] for r in current])):
            r['q' if key=='p' else key.replace('_p','_q')]=q

summary={assay:{'tested':sum(r['status']=='Tested' for r in results if r['assay']==assay),
                'nominal_p05':[r['gene'] for r in results if r['assay']==assay and r['p'] is not None and r['p']<.05],
                'fdr05':[r['gene'] for r in results if r['assay']==assay and r['q'] is not None and r['q']<.05]}
         for assay in ['rna','atac']}
notes=[
    'Groups are provisional marker-order patterns, not breakpoint-validated inversions. Translocations are not tested separately.',
    'Automatic assignments require SIRT6, CREB3L3 and TLE6 with one locus each on the same scaffold/chromosome. Coordinates are oriented along SIRT6, with both other marker TSSs downstream. Explicit user assignments override this rule; the original annotation result is preserved.',
    'All 31 non-human, non-macaque, non-bowhead species are documented. Excluded and unclassified species remain in the tables but never enter group tests. User-supplied assignments are identified separately from annotation evidence.',
    'The 63 original locus labels are the candidate set. No new gene alias mapping or orthology validation was applied to the historical expression/activity measurements.',
    'One species median is one observation; biological samples are not treated as independent species.',
    'For each target gene and assay, take classified species with a measurement for that target. Rank within each species using only candidate genes measured in ALL these species. Rank 1 is highest; ties receive average ranks.',
    'The shared ranking panel is target-specific. Rank denominators are identical between groups within a test, but can differ between target genes or assays. Missing values are never replaced by zero.',
    'At least 3 measured species per group and at least 2 shared genes are required to test. Constant ranks return p=1.',
    'Two-sided Wilcoxon rank-sum (Mann–Whitney U) p-values use exact conditional label permutations, computed by subset-sum dynamic programming with tied midranks. Two-sided p=2*min(lower-tail, upper-tail), capped at 1.',
    'BH false-discovery correction is performed separately for the tested genes in each assay. Sensitivity methods each have their own BH family.',
    'Positive rank-biserial effect means the target tends to rank nearer the top in rearranged species. Negative means nearer the bottom. This is a relative rank, not absolute expression or activity.',
    'Original-rank and relative-rank sensitivity tests use all measured original candidate genes within each species, so missing gene composition differs. Relative rank=(original rank-1)/(number measured-1). Raw-value tests are exploratory and do not normalize cross-species signal or source RNA values.',
    'Gene activity is species-median integrated CPM·bp in TSS ±5 kb. Highest-median locus was selected before the ≥50% available-window rule in the source dataset. RNA uses species-median source liver log.tpm without another log transform.',
    'No body mass, gene length or MLS adjustment is used: the predictor is gene-order group, not MLS. No phylogenetic adjustment is included in the Wilcoxon tests.',
    'Species share ancestry, and a single ancestral rearrangement can be inherited by multiple species. These exploratory p-values do not establish a causal inversion effect or independent evolutionary replication.',
    'Source RNA identity caveat: African spiny mouse species field says Acomys cahirinus, while the common-name field says Acomys russatus. The original species assignment is retained.',
    'Statistical tables are saved analysis results. Rerun analyze.py with the recorded assignment file and output directory, then the workbook builder, after changing groups, measurements or ranking rules.'
]
if overrides:
    notes.insert(2, 'User assignments, 31 October 2026: Transcaucasian mole vole is rearranged; eastern chipmunk, eastern cottontail, eastern mole, fox squirrel, naked mole-rat and North American porcupine are unrearranged. Muskrat is excluded. These override unresolved local annotations and are not independent breakpoint validation.')
payload=dict(species=classification,results=results,measurements=measurements,panels=panels,summary=summary,notes=notes,
             assignments=overrides,assignments_file=str(args.assignments.resolve()) if args.assignments else None,
             source=str(SOURCE),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             source_date='2026-10-29',analysis_date='2026-10-31',group_counts=dict(Counter(r['group'] for r in classification)),
             references=['https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html',
                         'https://www.journals.uchicago.edu/doi/10.1086/284325'])
(OUT/'analysis.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')

def write_tsv(name, records):
    with (OUT/name).open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]),delimiter='\t')
        writer.writeheader()
        for row in records:writer.writerow({k:'NA' if v is None else v for k,v in row.items()})
write_tsv('gene_rank_tests.tsv',results)
write_tsv('gene_rank_measurements.tsv',measurements)
flat=[]
for r in classification:
    item={k:v for k,v in r.items() if not k.endswith('_loci')}
    for gene in MARKERS:
        item[gene+'_coordinates']='; '.join(f'{l["chrom"]}:{l["start"]}-{l["end"]} ({l["strand"]})' for l in r[gene+'_loci']) or 'NA'
    flat.append(item)
# Uniform columns across classified and unclassified rows.
flat=[{key:r.get(key) for key in dict.fromkeys(k for rr in flat for k in rr)} for r in flat]
write_tsv('species_gene_order_groups.tsv',flat)
print(json.dumps({'groups':payload['group_counts'],'summary':summary,
                  'SIRT6':[r for r in results if r['gene']=='SIRT6']},indent=2))

os.environ.setdefault('MPLCONFIGDIR',str(OUT/'build/mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
colors={'Rearranged':'#0072B2','Unrearranged':'#D59B00'}
fig, axes=plt.subplots(1,2,figsize=(10,5.2),layout='constrained')
for ax, assay,title in zip(axes,['rna','atac'],['Liver RNA expression','TSS ±5 kb ATAC activity']):
    result=next(r for r in results if r['gene']=='SIRT6' and r['assay']==assay)
    obs=[r for r in measurements if r['gene']=='SIRT6' and r['assay']==assay and r['used_in_test']]
    for j,group in enumerate(GROUPS):
        rows=sorted([r for r in obs if r['group']==group],key=lambda r:r['rank'])
        ranks=[r['rank'] for r in rows]
        ax.boxplot(ranks,positions=[j],widths=.4,patch_artist=True,showfliers=False,
                   boxprops={'facecolor':colors[group]+'25','edgecolor':colors[group]},
                   medianprops={'color':'#222222','linewidth':1.4},whiskerprops={'color':colors[group]},capprops={'color':colors[group]})
        offsets=np.random.default_rng(52+j).uniform(-.12,.12,len(rows))
        ax.scatter(j+offsets,ranks,c=colors[group],s=35,edgecolor='white',linewidth=.6,zorder=3)
    ax.set_xticks([0,1],[f'Rearranged\n(n={result["n_rearranged"]})',f'Unrearranged\n(n={result["n_unrearranged"]})'])
    ax.set_ylim(result['reference_genes']+.5,.5)
    ax.set_yticks(range(1,result['reference_genes']+1) if result['reference_genes']<=10 else [1]+list(range(5,result['reference_genes']+1,5)))
    ax.set_ylabel(f'SIRT6 rank among {result["reference_genes"]} shared genes\n1 = highest')
    ax.set_title(f'{title}\nExact p = {result["p"]:.4g}; BH q = {result["q"]:.4g}',pad=14)
    ax.set_axisbelow(True);ax.yaxis.grid(True,color='#E5E7EB',linewidth=.7)
    for spine in ax.spines.values():spine.set_color('#222222')
fig.suptitle('SIRT6 rank by gene-order group',fontsize=14)
for suffix in ['png','svg','pdf']:fig.savefig(OUT/f'SIRT6_rank_comparison.{suffix}',dpi=240,facecolor='white')
plt.close(fig)

# Shared visual scale for target-specific rank effect, not a mixture of rank denominators.
fig,axes=plt.subplots(1,2,figsize=(12,10),layout='constrained')
for ax,assay,title in zip(axes,['rna','atac'],['Liver RNA expression','TSS ±5 kb ATAC activity']):
    rows=sorted([r for r in results if r['assay']==assay and r['status']=='Tested'],key=lambda r:r['delta_higher'])
    y=np.arange(len(rows));effect=[r['delta_higher'] for r in rows]
    ax.barh(y,effect,color=['#0072B2' if v>0 else '#D59B00' for v in effect],height=.75)
    ax.set_yticks(y,[r['gene']+(' **' if r['q']<.05 else ' *' if r['p']<.05 else '') for r in rows],fontsize=8)
    for tick,r in zip(ax.get_yticklabels(),rows):
        if r['gene']=='SIRT6':tick.set_fontweight('bold');tick.set_color('#A51624')
    ax.axvline(0,color='#333333',linewidth=.8);ax.set_xlim(-1.05,1.05)
    ax.set_title(title);ax.set_xlabel('Rank-biserial effect\nPositive = nearer the top in rearranged species')
    ax.set_axisbelow(True);ax.xaxis.grid(True,color='#E5E7EB',linewidth=.7)
fig.suptitle('Gene-rank associations with the rearranged pattern\n* nominal p < 0.05; ** BH q < 0.05. No phylogenetic adjustment.',fontsize=12)
for suffix in ['png','svg','pdf']:fig.savefig(OUT/f'Locus_gene_rank_effects.{suffix}',dpi=220,facecolor='white')
plt.close(fig)
