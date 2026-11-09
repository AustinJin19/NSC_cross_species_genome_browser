"""Historical positive ATAC-signal coverage audit, no longer used by live fits.

The denominator is the full requested region, including any span beyond the
scaffold. Only finite bigWig values > 0 count; zeros, missing bases and negative
values do not. Species coverage is the median across measured samples at the
selected locus. Scores and highest-activity locus selection are unchanged.
"""
import math
import statistics
from functools import lru_cache
from pathlib import Path

import pyBigWig

MIN_POSITIVE_SIGNAL_FRACTION = 0.5


@lru_cache(maxsize=8192)
def positive_bases(path, size, modified, chrom, start, end):
    """Count exact positive interval overlaps without expanding to base arrays."""
    from server import resolve_chrom
    with pyBigWig.open(path) as bw:
        sizes = bw.chroms()
        resolved = resolve_chrom(chrom, sizes)
        if resolved is None or end <= start:
            return None
        lo, hi = max(0, start), min(end, sizes[resolved])
        if hi <= lo:
            return 0
        return sum(min(hi, b) - max(lo, a)
                   for a, b, value in (bw.intervals(resolved, lo, hi) or ())
                   if math.isfinite(value) and value > 0)


def annotate_signal(row, chrom, start, end, folder):
    """Add per-sample and median species coverage; retain excluded measurements."""
    details, measured = {}, []
    requested = end - start
    for sample, value in row['samples'].items():
        info = dict(row.get('window_details', {}).get(sample, {}))
        path = Path(folder) / (sample + '.bw')
        positive = None
        if path.exists() and requested > 0:
            stat = path.stat()
            positive = positive_bases(str(path), stat.st_size, stat.st_mtime_ns,
                                      chrom, start, end)
        fraction = positive / requested if positive is not None else None
        info.update(positive_signal_bp=positive, positive_signal_fraction=fraction,
                    signal_eligible=fraction is not None and fraction >= MIN_POSITIVE_SIGNAL_FRACTION)
        details[sample] = info
        if value is not None and math.isfinite(value):
            measured.append(info)
    known = bool(measured) and all(d['positive_signal_fraction'] is not None for d in measured)
    fraction = statistics.median(d['positive_signal_fraction'] for d in measured) if known else None
    row.update(window_details=details,
               positive_signal_bp=statistics.median(d['positive_signal_bp'] for d in measured) if known else None,
               positive_signal_fraction=fraction,
               signal_eligible=fraction is not None and fraction >= MIN_POSITIVE_SIGNAL_FRACTION,
               signal_coverage_aggregation='median')
    return row
