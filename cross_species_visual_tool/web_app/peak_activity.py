"""Integrate sample bigWig signal over the union of its peaks inside a TSS window."""
import statistics
from functools import lru_cache
from pathlib import Path

import numpy as np
import pyBigWig

from regulatory_strength import peak_arrays


@lru_cache(maxsize=128)
def merged_peaks(path, size, modified):
    result = {}
    for chrom, (starts, ends) in peak_arrays(path, size, modified).items():
        # Sorted start/end events have zero active intervals between these gaps.
        # This also merges nested, duplicate and touching intervals.
        gaps = np.flatnonzero(starts[1:] > ends[:-1])
        result[chrom] = (starts[np.r_[0, gaps + 1]], ends[np.r_[gaps, len(ends) - 1]])
    return result


@lru_cache(maxsize=8192)
def window_score(bw_path, bw_size, bw_modified, bed_path, bed_size, bed_modified,
                 chrom, start, end):
    from server import resolve_chrom
    from calculate_gene_activity import integrate

    details = dict(status='chromosome_unavailable', quantified_bp=0, window_bp=0,
                   peak_bp=0, covered_bp=0)
    with pyBigWig.open(bw_path) as bw:
        sizes = bw.chroms()
        c = resolve_chrom(chrom, sizes)
        if c is None:
            return None, details
        lo, hi = max(0, start), min(sizes[c], end)
        if hi <= lo:
            details['status'] = 'window_outside_chromosome'
            return None, details
        details['window_bp'] = hi - lo
        groups = merged_peaks(bed_path, bed_size, bed_modified)
        pc = resolve_chrom(chrom, groups)
        if pc is None and groups and not any(resolve_chrom(k, sizes) for k in groups):
            details['status'] = 'peak_assembly_mismatch'
            return None, details
        details['status'] = 'clipped' if lo != start or hi != end else 'ok'
        if pc is None:
            return 0., details
        starts, ends = groups[pc]
        first = np.searchsorted(ends, lo, side='right')
        last = np.searchsorted(starts, hi, side='left')
        if first >= last:
            return 0., details
        starts = np.maximum(starts[first:last], lo)
        ends = np.minimum(ends[first:last], hi)
        details['quantified_bp'] = details['peak_bp'] = int(np.sum(ends - starts))
        values, coverage = integrate(bw.intervals(c, int(starts[0]), int(ends[-1])), starts, ends)
        details['covered_bp'] = int(np.sum(coverage))
        return float(np.sum(values)), details


def score_records(records, flank):
    from server import SPECS, peak_path
    root = Path(__file__).resolve().parent.parent
    folders = {s[0]: root / s[4] for s in SPECS}
    result = []
    for original in records:
        row = dict(original)
        chrom, span, strand = row['locus'].rsplit('|', 1)[1].rsplit(':', 2)
        a, b = map(int, span.split('-'))
        tss = a if strand == '+' else b - 1
        samples, details = {}, {}
        for sample in row['samples']:
            bw = folders[row['species']] / (sample + '.bw')
            bed = peak_path(bw)
            if not bw.exists() or not bed.exists():
                samples[sample] = None
                details[sample] = dict(status='signal_file_missing' if not bw.exists() else 'peak_file_missing',
                                       quantified_bp=0, window_bp=0, peak_bp=0, covered_bp=0)
                continue
            bs, ps = bw.stat(), bed.stat()
            samples[sample], details[sample] = window_score(
                str(bw), bs.st_size, bs.st_mtime_ns, str(bed), ps.st_size, ps.st_mtime_ns,
                chrom, tss - flank, tss + flank)
        valid = [v for v in samples.values() if v is not None]
        # Preserve source units: human tracks are scaled signal, not CPM.
        row.update(samples=samples, median=statistics.median(valid) if valid else None,
                   window_details=details, region_start=tss-flank, region_end=tss+flank)
        result.append(row)
    return result
