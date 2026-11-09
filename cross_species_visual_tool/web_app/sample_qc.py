"""Shared, portable sample exclusions for ATAC tracks and activity calculations."""
import json
import math
import statistics
from pathlib import Path

QC_PATH = Path(__file__).resolve().parent / 'data/sample_qc_exclusions.json'
QC_POLICY = json.loads(QC_PATH.read_text())


def exclusions_for_policy(policy):
    """Allow the TSS filter to be disabled without changing legacy exclusions."""
    excluded = set(policy.get('always_exclude', []))
    if policy.get('enabled', True):
        excluded.update(row['sample'] for row in policy['excluded_samples'])
    return frozenset(excluded)


EXCLUDED_SAMPLES = exclusions_for_policy(QC_POLICY)


def sample_id(sample):
    """Normalize browser IDs and indexed bigWig stems to the same QC sample ID."""
    return Path(str(sample)).name.removesuffix('.bw').removesuffix('.cpm').removesuffix('_scaled')


def filter_records(records):
    """Filter before scoring/choosing loci; rebuild stored body medians as well.

    Keep unavailable measurements as None. A species with no retained samples
    contributes no records, never a zero measurement. Original records are kept
    intact for reproducibility.
    """
    result = []
    for original in records:
        row = dict(original)
        samples = json.loads(row['samples']) if isinstance(row['samples'], str) else row['samples']
        row['samples'] = {name: value for name, value in samples.items()
                          if sample_id(name) not in EXCLUDED_SAMPLES}
        if not row['samples']:
            continue
        valid = [v for v in row['samples'].values() if v is not None and math.isfinite(v)]
        row['median'] = statistics.median(valid) if valid else None
        if 'window_details' in row:
            row['window_details'] = {name: detail for name, detail in row['window_details'].items()
                                     if name in row['samples']}
        result.append(row)
    return result
