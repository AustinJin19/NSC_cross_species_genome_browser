"""Adult liver dbSUPER calls on native hg38/mm10, exposed as 1-based intervals."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from functools import lru_cache
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[2] / 'data/regulatory_annotations/dbSUPER_liver'
SPECIES = {'Human': ('human_hg38', 'hg38'), 'Mouse': ('mouse_mm10', 'mm10')}


def version(path):
    return path.stat().st_mtime_ns if path.is_file() else None


def bed_rows(path):
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            chrom, start, end, name = line.rstrip('\n').split('\t')[:4]
            # The browser uses 1-based inclusive coordinates; source BED is 0-based.
            yield chrom, int(start) + 1, int(end), name


@lru_cache(maxsize=4)
def read_index(directory, sid, all_version, candidate_version, metadata_version):
    directory = Path(directory)
    prefix, _ = SPECIES[sid]
    candidates = set(bed_rows(directory / f'{prefix}_candidate_liver_specific_super_enhancers.bed')) if candidate_version is not None else set()
    metadata = json.loads((directory / 'metadata.json').read_text()) if metadata_version is not None else {}
    species_metadata = metadata.get('species', {}).get(sid, {})
    region_metadata = species_metadata.get('regions', {})
    chroms = defaultdict(list)
    rows = list(bed_rows(directory / f'{prefix}_liver_super_enhancers.bed'))
    if not candidates.issubset(set(rows)):
        raise ValueError('Candidate super-enhancer BED must be a subset of the liver BED')
    for chrom, start, end, name in rows:
        if start < 1 or end < start:
            raise ValueError('Invalid super-enhancer interval')
        info = region_metadata.get(name, {})
        chroms[chrom].append(dict(start=start, end=end, name=name,
                                 candidateLiverSpecific=(chrom, start, end, name) in candidates,
                                 nLiverDatasets=info.get('nLiverDatasets'),
                                 liverDatasets=info.get('liverDatasets', []),
                                 sourceCalls=info.get('sourceCalls', [])))
    index = {}
    for chrom, entries in chroms.items():
        entries.sort(key=lambda r: (r['start'], r['end']))
        ends, stop = [], 0
        for row in entries:
            stop = max(stop, row['end'])
            ends.append(stop)
        index[chrom] = ([r['start'] for r in entries], ends, entries)
    return index, len(rows), len(candidates), species_metadata


def query(sid, chrom, lo, hi):
    if sid not in SPECIES:
        return None
    prefix, assembly = SPECIES[sid]
    all_version = version(DATA / f'{prefix}_liver_super_enhancers.bed')
    candidate_version = version(DATA / f'{prefix}_candidate_liver_specific_super_enhancers.bed')
    result = dict(available=all_version is not None, candidateAvailable=candidate_version is not None,
                  assembly=assembly, source='dbSUPER v2.0', sourceUrl='https://dbsuper.khanlab.bio/',
                  snapshot='2026-11-02', regions=[], total=0, candidateTotal=0,
                  color='#7B3294', datasets=[])
    if not result['available']:
        return result
    index, total, candidate_total, metadata = read_index(
        str(DATA), sid, all_version, candidate_version, version(DATA / 'metadata.json'))
    result.update(total=total, candidateTotal=candidate_total, datasets=metadata.get('datasets', []),
                  comparatorDatasets=metadata.get('summary', {}).get('adult_nonliver_datasets'))
    if chrom not in index or hi < max(1, lo):
        return result
    starts, ends, rows = index[chrom]
    result['regions'] = [r for r in rows[bisect_left(ends, lo):bisect_right(starts, hi)] if r['end'] >= lo]
    return result
