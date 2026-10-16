"""Count BED peaks overlapping half-open TSS windows; no distance weighting."""
import gzip,statistics
from pathlib import Path
from functools import lru_cache
import numpy as np
import pyBigWig

@lru_cache(maxsize=128)
def peak_arrays(path, size, modified):
    groups={}
    with (gzip.open(path,'rt') if path.endswith('.gz') else open(path)) as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            c,a,b,*_=line.split();a,b=int(a),int(b)
            if b>a:groups.setdefault(c,[]).append((a,b))
    return {c:(np.sort([a for a,b in rows]),np.sort([b for a,b in rows])) for c,rows in groups.items()}

@lru_cache(maxsize=128)
def chrom_sizes(path,size,modified):
    with pyBigWig.open(path) as bw:return bw.chroms()

def score_records(records,flank):
    from server import SPECS,resolve_chrom,peak_path
    root=Path(__file__).resolve().parent.parent
    folders={s[0]:root/s[4] for s in SPECS}
    output=[]
    for original in records:
        row=dict(original);chrom,span,strand=row['locus'].rsplit('|',1)[1].rsplit(':',2)
        a,b=map(int,span.split('-'));tss=a if strand=='+' else b-1
        samples,details={},{}
        for sample in row['samples']:
            bw=folders[row['species']]/(sample+'.bw')
            stem=sample.removesuffix('.cpm').removesuffix('_scaled')
            path=peak_path(bw)
            value=None;width=0;status='peak_file_missing'
            if path.exists() and bw.exists():
                stat=bw.stat();sizes=chrom_sizes(str(bw),stat.st_size,stat.st_mtime_ns);c=resolve_chrom(chrom,sizes)
                if c is None:status='chromosome_unavailable'
                else:
                    lo=max(0,tss-flank);hi=min(sizes[c],tss+flank);width=max(0,hi-lo)
                    stat=path.stat();groups=peak_arrays(str(path),stat.st_size,stat.st_mtime_ns);pc=resolve_chrom(chrom,groups)
                    # No peaks on a supported chromosome is a real zero.
                    if pc is None and groups and not any(resolve_chrom(k,sizes) for k in groups):status='peak_assembly_mismatch'
                    elif width<=0:status='window_outside_chromosome'
                    else:
                        value=0 if pc is None else int(np.searchsorted(groups[pc][0],hi,side='left')-np.searchsorted(groups[pc][1],lo,side='right'))
                        status='clipped' if width<2*flank else 'ok'
            samples[sample]=value;details[sample]={'status':status,'quantified_bp':width}
        valid=[v for v in samples.values() if v is not None]
        row.update(samples=samples,median=statistics.median(valid) if valid else None,unit='peak_count',window_details=details,region_start=tss-flank,region_end=tss+flank)
        output.append(row)
    return output
