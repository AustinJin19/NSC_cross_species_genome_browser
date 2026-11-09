"""Species-level ATAC candidate SE shared (>=2 samples) and union (>=1) regions.

Installed BED files are 0-based, half-open; query coordinates and responses use
the browser's 1-based, inclusive convention. These are independent of dbSUPER.
"""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from functools import lru_cache
import json
from pathlib import Path

from super_enhancers import bed_rows, version

DATA = Path(__file__).resolve().parents[2] / 'data/regulatory_annotations/atac_super_enhancers'


@lru_cache(maxsize=4)
def metadata(directory, stamp):
    return json.loads((Path(directory) / 'metadata.json').read_text()) if stamp is not None else {}


@lru_cache(maxsize=64)
def read_index(path, stamp):
    chroms = defaultdict(list)
    for chrom, start, end, name in bed_rows(Path(path)):
        if start < 1 or end < start:
            raise ValueError('Invalid ATAC candidate super-enhancer interval')
        chroms[chrom].append(dict(start=start, end=end, name=name))
    index = {}
    for chrom, rows in chroms.items():
        rows.sort(key=lambda r: (r['start'], r['end']))
        ends, stop = [], 0
        for row in rows:
            stop = max(stop, row['end'])
            ends.append(stop)
        index[chrom] = ([r['start'] for r in rows], ends, rows)
    return index


def query(sid, chrom, lo, hi, region_set='shared_ge2'):
    if region_set not in ('shared_ge2', 'union'):
        raise ValueError('Unknown ATAC super-enhancer region set')
    info = metadata(str(DATA), version(DATA / 'metadata.json')).get('species', {}).get(sid)
    result = dict(available=False, status='not_provided', source='ATAC candidate super-enhancers',
                  regionSet=region_set, minimumSupport=1 if region_set == 'union' else 2,
                  regions=[], total=0, nSamples=None, samples=[], color='#B65C00')
    if info is None:
        return result
    result.update(nSamples=info['nSamples'], samples=info['samples'])
    selected = info.get('union') if region_set == 'union' else info
    if selected is None:
        return result
    path = DATA / selected['file']
    stamp = version(path)
    result.update(available=stamp is not None,
                  status=selected['supportStatus'] if stamp is not None else 'file_unavailable',
                  total=selected['regions'], file=selected['file'],
                  coveredBp=selected.get('coveredBp', selected.get('sharedBp')))
    if region_set == 'shared_ge2':
        result['sharedBp'] = selected['sharedBp']
    if stamp is None or hi < max(1, lo):
        return result
    index = read_index(str(path), stamp)
    if chrom in index:
        starts, ends, rows = index[chrom]
        result['regions'] = [r for r in rows[bisect_left(ends, lo):bisect_right(starts, hi)] if r['end'] >= lo]
    return result
