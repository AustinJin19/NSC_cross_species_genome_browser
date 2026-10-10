"""Hypothalamus ABC predictions, independent of the selected liver replicate."""
from pathlib import Path
from functools import lru_cache
import gzip
import math

SOURCE = Path(__file__).resolve().parents[2] / 'data/interactions/nmr_abc/mHetGlaV3.primary_hypothalamus_ABC_enhancer.bedpe (1).gz'

@lru_cache(maxsize=1)
def read_links():
    if not SOURCE.exists():
        return None
    rows = []
    with gzip.open(SOURCE, 'rt') as handle:
        next(handle)  # Named header, not a BED record.
        for line in handle:
            f = line.rstrip().split('\t')
            a,b,c,d = map(int,(f[1],f[2],f[4],f[5]))
            score = float(f[7])
            if a < 0 or b <= a or c < 0 or d < c or not math.isfinite(score):
                raise ValueError('Invalid ABC record')
            # Source target anchors have start=end. Treat as a BED point,
            # extending to one base for display; preserve original coordinates.
            rows.append(dict(chrom1=f[0], start1=a+1, end1=b, chrom2=f[3],
                start2=c+1, end2=max(c+1,d), name=f[6], score=score,
                targetGene=f[6].split('_',1)[0], sourceTargetStart=c,sourceTargetEnd=d))
    return rows


def query(chrom, lo, hi, limit=300):
    rows = read_links()
    hits = [r for r in (rows or []) if
        (r['chrom1']==chrom and r['start1']<=hi and r['end1']>=lo) or
        (r['chrom2']==chrom and r['start2']<=hi and r['end2']>=lo)]
    hits.sort(key=lambda r:(-r['score'],r['name']))
    return dict(available=rows is not None, regions=hits[:limit],total=len(hits),
                label='ABC predictions · hypothalamus · TEST',kind='abc')
