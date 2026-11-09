"""Genomic-span eligibility shared by activity and lifespan analyses.

Availability is the requested interval intersected with the matching bigWig
chromosome/scaffold bounds, not read coverage or the number of bases in peaks.
Raw measurements are retained for inspection; fits must check region_eligible
after choosing the highest-activity locus for each species.
"""
import math
from functools import lru_cache
from pathlib import Path

import pyBigWig

MIN_AVAILABLE_FRACTION = 0.5


@lru_cache(maxsize=128)
def chrom_sizes(path, size, modified):
    with pyBigWig.open(path) as bw:
        return bw.chroms()


def coverage(chrom, start, end, sizes):
    from server import resolve_chrom
    requested = end - start
    resolved = resolve_chrom(chrom, sizes)
    available = (max(0, min(end, sizes[resolved]) - max(0, start))
                 if resolved is not None and requested > 0 else 0)
    fraction = available / requested if requested > 0 else None
    return dict(requested_bp=requested, available_bp=available,
                available_fraction=fraction,
                region_eligible=requested > 0 and available * 2 >= requested)


def body_region(locus):
    chrom, span, strand = locus.rsplit('|', 1)[1].rsplit(':', 2)
    start, end = map(int, span.split('-'))
    return (chrom, start - 2000, end) if strand == '+' else (chrom, start, end + 2000)


def annotate_region(row, chrom, start, end, folder):
    """Attach coverage without changing scores or multiple-locus selection.

    All measured samples must meet the threshold. Normally their chromosome
    sizes are identical; using the minimum is conservative if they disagree.
    Unknown bounds for a measured sample cannot establish eligibility.
    """
    details, measured = {}, []
    for sample, value in row['samples'].items():
        info = dict(row.get('window_details', {}).get(sample, {}))
        path = Path(folder) / (sample + '.bw')
        if path.exists():
            stat = path.stat()
            info.update(coverage(chrom, start, end,
                                 chrom_sizes(str(path), stat.st_size, stat.st_mtime_ns)))
        else:
            info.update(requested_bp=end-start, available_bp=None,
                        available_fraction=None, region_eligible=False)
        details[sample] = info
        if value is not None and math.isfinite(value):
            measured.append(info)
    known = measured and all(d['available_fraction'] is not None for d in measured)
    limiting = min(measured, key=lambda d: d['available_fraction']) if known else None
    row.update(window_details=details, region_start=start, region_end=end,
               requested_bp=end-start,
               available_bp=limiting['available_bp'] if limiting else None,
               available_fraction=limiting['available_fraction'] if limiting else None,
               region_eligible=bool(measured) and all(d['region_eligible'] for d in measured))
    return row
