"""Liver ChIP-seq BED regions, converted to browser 1-based inclusive coordinates."""
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
from bisect import bisect_left,bisect_right
ROOT=Path(__file__).resolve().parents[2]/'mhetglav3'
TRACKS=[('H3K4me3','#d62728','hetGla-H3K4me3_replicated-peaks_mHetGlaV3.bed'),('H3K27ac','#e6b800','hetGla-H3K27Ac_replicated-peaks_mHetGlaV3.bed')]

LIFTED=ROOT.parent/'ChIP_seq'/'liftover'/'final'
ASSEMBLIES={'Human':'hg38','Mouse':'mm10','Rat':'GRCr8','Macaque':'T2T-MMU8v2.0'}

@lru_cache(maxsize=16)
def index(path,version):
    groups=defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#','track ','browser ')):continue
            fields=line.split();start,end=int(fields[1]),int(fields[2])
            if start<0 or end<=start:raise ValueError(f'Invalid ChIP interval in {path.name}')
            groups[fields[0]].append((start+1,end,fields[3] if len(fields)>3 else '.'))
    out={}
    for chrom,rows in groups.items():
        rows.sort();maximum=0;prefix=[]
        for start,end,_ in rows:maximum=max(maximum,end);prefix.append(maximum)
        out[chrom]=(rows,[r[0] for r in rows],prefix)
    return out

def query(chrom,lo,hi,species='NMR'):
    tracks=[]
    if species!='NMR' and species not in ASSEMBLIES:return tracks
    assembly='mHetGlaV3' if species=='NMR' else ASSEMBLIES[species]
    for name,color,filename in TRACKS:
        if species=='NMR':path=ROOT/filename
        else:
            mark='H3K27Ac' if name=='H3K27ac' else name
            filename=f'{species}_{mark}_{assembly}.bed'
            path=LIFTED/filename
        regions=[]
        if path.exists():
            data=index(path,path.stat().st_mtime_ns)
            if chrom in data:
                rows,starts,prefix=data[chrom]
                regions=[dict(start=a,end=b,name=n) for a,b,n in rows[bisect_left(prefix,lo):bisect_right(starts,hi)] if b>=lo]
        tracks.append(dict(name=name,color=color,source=filename,assembly=assembly,available=path.exists(),regions=regions))
    return tracks
