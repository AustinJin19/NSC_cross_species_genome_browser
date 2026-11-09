"""Indexed human liver DeepLoop contacts, lifted to hg38 BEDPE."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from functools import lru_cache
import gzip
from pathlib import Path

LOOPS = Path(__file__).resolve().parents[2] / 'data/interactions/human_liver_deeploop/GSE167200_Liver.top300K.hg38.bedpe.gz'


@lru_cache(maxsize=2)
def read_index(path, version):
    """Index both anchors; prefix maximum ends also find long overlapping bins."""
    rows, chroms = [], defaultdict(list)
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            row = (f[0], int(f[1])+1, int(f[2]), f[3], int(f[4])+1, int(f[5]), f[6], float(f[7]), f[10], f[11])
            index = len(rows)
            rows.append(row)
            chroms[row[0]].append((row[1], row[2], index))
            chroms[row[3]].append((row[4], row[5], index))
    result = {}
    for chrom, entries in chroms.items():
        entries.sort()
        ends, end = [], 0
        for _, stop, _ in entries:
            end = max(end, stop)
            ends.append(end)
        result[chrom] = ([entry[0] for entry in entries], ends, entries)
    return rows, result


def query(chrom, lo, hi, limit=300):
    result = dict(available=False, regions=[], total=0, limit=limit, kind='deeploop',
                  label='Liver Hi-C loops · DeepLoop', assembly='hg38',
                  source='GSE167200 · liver', sourceUrl='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167200')
    if not LOOPS.is_file():
        return result
    rows, chroms = read_index(str(LOOPS), LOOPS.stat().st_mtime_ns)
    result['available'] = True
    if chrom not in chroms or hi < lo:
        return result
    starts, ends, entries = chroms[chrom]
    hits = {index for start, end, index in entries[bisect_left(ends, lo):bisect_right(starts, hi)] if end >= lo}
    ordered = sorted(hits, key=lambda index: (-rows[index][7], rows[index][6]))
    keys = ('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'hg19Anchor1', 'hg19Anchor2')
    result.update(total=len(ordered), regions=[dict(zip(keys, rows[index])) for index in ordered[:limit]])
    return result
