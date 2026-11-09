"""Chromosome/length-preserving permutation test of unique LAR participation."""
from pathlib import Path
from collections import defaultdict
import argparse,csv,json,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
METRICS=['LAR_any_interaction','LAR_one_hit_interaction','LAR_two_hit_interaction','LAR_only_one_hit','LAR_only_two_hit','LAR_both_types','LAR_no_interaction','interactions_one_hit','interactions_two_hit']

def load():
    sizes={r[0]:int(r[1]) for r in (l.split() for l in (ROOT/'reference_genomes/hg19/hg19.chrom.sizes').open())}
    lars=[l.split() for l in (ROOT/'data/regulatory_annotations/pgls/pgls_LARs_hg19.bed').open() if l.strip() and not l.startswith('#')]
    assert len({r[3] for r in lars})==len(lars)
    for c,a,b,*_ in lars:assert c in sizes and 0<=int(a)<int(b)<=sizes[c],(c,a,b)
    anchors=[];lookup={};edges=[]
    with (ROOT/'data/interactions/hicorr/Merged.raw_HiCorr_DeepLoop.1.5.bedpe').open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            r=line.split();pair=[]
            for i in (0,3):
                key=(r[i],int(r[i+1]),int(r[i+2]))
                assert key[0] in sizes and 0<=key[1]<key[2]<=sizes[key[0]],key
                if key not in lookup:lookup[key]=len(anchors);anchors.append(key)
                pair.append(lookup[key])
            edges.append(pair)
    bychrom=defaultdict(list)
    for i,(c,a,b) in enumerate(anchors):bychrom[c].append((a,b,i))
    index={}
    for c,rows in bychrom.items():
        rows=np.array(sorted(rows),dtype=np.int64);index[c]=(rows[:,0],rows[:,1],np.maximum.accumulate(rows[:,1]),rows[:,2])
    groups={c:np.array([i for i,r in enumerate(lars) if r[0]==c]) for c in {r[0] for r in lars}}
    return sizes,lars,anchors,np.array(edges,dtype=np.int32),index,groups

def pairs_for_regions(starts,lengths,index,groups):
    all_l=[];all_a=[]
    for c,ids in groups.items():
        if c not in index:continue
        a,b,prefix,anchorids=index[c];s=starts[ids];e=s+lengths[ids]
        lo=np.searchsorted(prefix,s,side='right');hi=np.searchsorted(a,e,side='left');counts=hi-lo
        keep=counts>0
        ids=ids[keep];lo=lo[keep];counts=counts[keep]
        if not len(ids):continue
        offsets=np.cumsum(counts)-counts
        q=np.repeat(lo,counts)+np.arange(counts.sum())-np.repeat(offsets,counts)
        larids=np.repeat(ids,counts);valid=b[q]>starts[larids]
        all_l.append(larids[valid]);all_a.append(anchorids[q[valid]])
    return (np.concatenate(all_l),np.concatenate(all_a)) if all_l else (np.array([],dtype=int),np.array([],dtype=int))

def score(starts,lengths,index,groups,edges,nanchors,details=False):
    lids,aids=pairs_for_regions(starts,lengths,index,groups)
    active=np.zeros(nanchors,dtype=np.uint8);active[aids]=1
    edge_hits=active[edges[:,0]]+active[edges[:,1]]
    flags=np.zeros(nanchors,dtype=np.uint8)
    one=edges[edge_hits==1].ravel();two=edges[edge_hits==2].ravel()
    flags[one]|=1;flags[two]|=2
    larflags=np.zeros(len(starts),dtype=np.uint8)
    np.bitwise_or.at(larflags,lids,flags[aids])
    counts=[np.count_nonzero(larflags),np.count_nonzero(larflags&1),np.count_nonzero(larflags&2),np.count_nonzero(larflags==1),np.count_nonzero(larflags==2),np.count_nonzero(larflags==3),np.count_nonzero(larflags==0),np.count_nonzero(edge_hits==1),np.count_nonzero(edge_hits==2)]
    return (np.array(counts),larflags) if details else np.array(counts)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--permutations',type=int,default=1000);ap.add_argument('--seed',type=int,default=20260918);args=ap.parse_args()
    started=time.time();sizes,lars,anchors,edges,index,groups=load()
    starts=np.array([int(r[1]) for r in lars]);lengths=np.array([int(r[2])-int(r[1]) for r in lars]);limits=np.array([sizes[r[0]] for r in lars])-lengths+1
    observed,flags=score(starts,lengths,index,groups,edges,len(anchors),True)
    # Independent verification against the previously generated annotated BEDPE.
    one=set();two=set()
    with (HERE.parent/'HiCorr_pgls_LARs_hg19.bedpe').open() as f:
        next(f)
        for line in f:
            r=line.rstrip().split('\t');names=set((r[10]+';'+r[14]).split(';'))-{'.'}
            (two if r[8]=='both' else one).update(names)
    assert {r[3] for r,v in zip(lars,flags) if v&1}==one
    assert {r[3] for r,v in zip(lars,flags) if v&2}==two
    print('Observed',dict(zip(METRICS,map(int,observed))),flush=True)
    with (HERE/'LAR_interaction_membership.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['chrom','start','end','LAR_ID','any_interaction','one_hit_interaction','two_hit_interaction','category'])
        for row,v in zip(lars,flags):w.writerow(row[:4]+[int(v>0),int(bool(v&1)),int(bool(v&2)),['none','one_hit_only','two_hit_only','both_types'][v]])
    rng=np.random.default_rng(args.seed);null=[]
    with (HERE/'shuffle_counts.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['permutation']+METRICS)
        for i in range(args.permutations):
            shuffled=rng.integers(0,limits)
            values=score(shuffled,lengths,index,groups,edges,len(anchors));null.append(values);w.writerow([i+1]+values.tolist())
            if (i+1)%100==0:print(f'{i+1}/{args.permutations} permutations; {time.time()-started:.1f}s',flush=True)
    null=np.array(null);summary=[]
    for j,m in enumerate(METRICS):
        v=null[:,j];obs=int(observed[j]);mean=float(v.mean())
        summary.append(dict(metric=m,observed=obs,observed_percent_of_LARs=100*obs/len(lars) if m.startswith('LAR_') else None,null_mean=mean,null_sd=float(v.std(ddof=1)),null_2_5_percentile=float(np.quantile(v,.025)),null_97_5_percentile=float(np.quantile(v,.975)),fold_enrichment=obs/mean if mean else None,permutations_at_least_observed=int(np.sum(v>=obs)),p_enrichment=float((1+np.sum(v>=obs))/(args.permutations+1)),p_depletion=float((1+np.sum(v<=obs))/(args.permutations+1))))
    with (HERE/'enrichment_summary.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
    metadata=dict(permutations=args.permutations,seed=args.seed,n_LARs=len(lars),n_interactions=len(edges),n_unique_anchors=len(anchors),LARs_on_contigs_without_HiCorr=sum(len(ids) for c,ids in groups.items() if c not in index),chrom_sizes_source='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',null='Independently sample each start uniformly from 0 through chromosome_length-region_length inclusive; preserve original chromosome and exact length; allow overlap between shuffled regions; no gap/blacklist/accessibility mask',coordinates='hg19 BED 0-based half-open; >=1 bp overlap',elapsed_seconds=time.time()-started)
    metadata['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'data/regulatory_annotations/pgls/pgls_LARs_hg19.bed',ROOT/'data/interactions/hicorr/Merged.raw_HiCorr_DeepLoop.1.5.bedpe',ROOT/'reference_genomes/hg19/hg19.chrom.sizes']}
    (HERE/'analysis_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
