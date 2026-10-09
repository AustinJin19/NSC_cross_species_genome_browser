"""Sample-specific BEDPE loop retrieval; coordinates returned 1-based inclusive."""
from pathlib import Path
from functools import lru_cache
import math

LOOPS = Path(__file__).resolve().parents[2] / 'loops'

@lru_cache(maxsize=3)
def read_loops(sample):
    if sample not in {'NMR3Liver', 'NMR4Liver', 'NMR8Liver'}:
        return None
    path = LOOPS / f'{sample}.microc_ch0.enrich1.0.loops.bedpe'
    if not path.exists():
        return None
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            a, b, c, d = map(int, (f[1], f[2], f[4], f[5]))
            if a < 0 or c < 0 or b <= a or d <= c:
                raise ValueError(f'Invalid loop coordinates in {path.name}')
            score, oe, enrichment = map(float, (f[7], f[10], f[11]))
            if not all(math.isfinite(v) for v in (score, oe, enrichment)):
                raise ValueError(f'Invalid loop score in {path.name}')
            rows.append(dict(chrom1=f[0], start1=a+1, end1=b, chrom2=f[3],
                             start2=c+1, end2=d, name=f[6], score=score,
                             log2OE=oe, log2Enrichment=enrichment,
                             windowsCalled=int(f[12]), windowsCovering=int(f[13])))
    return rows


def query(sample, chrom, lo, hi, limit=300):
    rows = read_loops(sample)
    if rows is None:
        return dict(available=False, regions=[], total=0, sample=sample)
    # Include a loop if either anchor overlaps the window. Retain the remote
    # anchor's original coordinates for accurate off-screen labels.
    hits = [r for r in rows if
            (r['chrom1'] == chrom and r['start1'] <= hi and r['end1'] >= lo) or
            (r['chrom2'] == chrom and r['start2'] <= hi and r['end2'] >= lo)]
    hits.sort(key=lambda r: (-r['score'], r['name']))
    return dict(available=True, regions=hits[:limit], total=len(hits), sample=sample,
                limit=limit, label='Predicted Micro-C loops · enrich1.0')
