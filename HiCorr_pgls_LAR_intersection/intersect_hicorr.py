"""Annotate either BEDPE anchor with overlapping hg19 pgls LARs (BED half-open)."""
from pathlib import Path
from collections import defaultdict, Counter
from bisect import bisect_left, bisect_right
from functools import lru_cache
import json

ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
LARS=ROOT/'data/regulatory_annotations/pgls/pgls_LARs_hg19.bed'
LOOPS=ROOT/'data/interactions/hicorr/Merged.raw_HiCorr_DeepLoop.1.5.bedpe'

def index_lars(path):
    intervals=defaultdict(list)
    with path.open() as f:
        for number,line in enumerate(f,1):
            if not line.strip() or line.startswith(('#','track ','browser ')): continue
            row=line.rstrip('\r\n').split('\t')
            chrom,start,end,name,score=row
            start,end=int(start),int(end)
            if start<0 or end<=start:raise ValueError(f'Invalid LAR interval at row {number}')
            intervals[chrom].append((start,end,name))
    result={}
    for chrom,rows in intervals.items():
        rows.sort();maximum=0;prefix=[]
        for start,end,name in rows:
            maximum=max(maximum,end);prefix.append(maximum)
        result[chrom]=(rows,[r[0] for r in rows],prefix)
    return result

def overlap(index,chrom,start,end):
    if chrom not in index:return ()
    rows,starts,prefix=index[chrom]
    lower=bisect_right(prefix,start);upper=bisect_left(starts,end)
    return tuple((name,f'{chrom}:{a}-{b}',min(end,b)-max(start,a)) for a,b,name in rows[lower:upper] if b>start)

def main():
    index=index_lars(LARS)
    @lru_cache(maxsize=400000)
    def hits(chrom,start,end):return overlap(index,chrom,start,end)
    def columns(values):
        return [str(len(values)),';'.join(v[0] for v in values) or '.', ';'.join(v[1] for v in values) or '.', ';'.join(str(v[2]) for v in values) or '.']
    counts=Counter();seen=set();anchor_seen=set()
    target=OUT/'HiCorr_pgls_LARs_hg19.bedpe'
    header=['chrom1','start1','end1','chrom2','start2','end2','HiCorr_ID','HiCorr_score','pgls_LAR_overlap']
    for side in (1,2):header += [f'anchor{side}_pgls_LAR_count',f'anchor{side}_pgls_LAR_IDs',f'anchor{side}_pgls_LAR_hg19_intervals',f'anchor{side}_overlap_bp']
    with LOOPS.open() as src,target.open('w') as dst:
        dst.write('#'+'\t'.join(header)+'\n')
        for number,line in enumerate(src,1):
            if not line.strip() or line.startswith(('#','track ','browser ')):continue
            row=line.rstrip('\r\n').split('\t')
            if len(row)!=8:raise ValueError(f'Expected 8 HiCorr columns at row {number}')
            c1,s1,e1,c2,s2,e2=row[:6];s1,e1,s2,e2=map(int,(s1,e1,s2,e2))
            if min(s1,s2)<0 or e1<=s1 or e2<=s2:raise ValueError(f'Invalid BEDPE coordinates at row {number}')
            a,b=hits(c1,s1,e1),hits(c2,s2,e2);counts['input_interactions']+=1
            if not a and not b:continue
            side='both' if a and b else 'anchor1' if a else 'anchor2'
            counts[side]+=1;counts['output_interactions']+=1
            for anchor,values in [((c1,s1,e1),a),((c2,s2,e2),b)]:
                if values:anchor_seen.add(anchor)
                seen.update(v[0] for v in values)
            dst.write('\t'.join(row+[side]+columns(a)+columns(b))+'\n')
    counts['input_pgls_LARs']=sum(len(v[0]) for v in index.values())
    counts['overlapped_unique_pgls_LARs']=len(seen)
    counts['overlapped_unique_HiCorr_anchors']=len(anchor_seen)
    summary=dict(counts,assembly='hg19 (both inputs; HiCorr confirmed by user)',coordinates='0-based, half-open',rule='Keep each HiCorr row if either anchor overlaps any pgls LAR by at least 1 bp',output=target.name)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
