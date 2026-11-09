#!/usr/bin/env python3
"""Count the canonical per-sample narrowPeak records; retain every QC level."""
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
import statistics
import sys

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
MANIFEST = ROOT / 'data/liver_atac_bundle/manifest.local.tsv'
EXCLUDE = {'Homo_sapiens', 'Macaca_mulatta', 'Balaena_mysticetus'}
read_tsv = lambda p: list(csv.DictReader(p.open(), delimiter='\t'))
common = {r['species']: r['common_name'] for r in read_tsv(ROOT/'tmm_pca/all_genes/sample_correlation/sample_order.tsv')}
qc = {r['sample']: r for r in read_tsv(ROOT/'bundle_sample_qc_summary.tsv')}
manifest = read_tsv(MANIFEST)
rows = [r for r in manifest if r['species'] not in EXCLUDE]
assert len({r['sample'] for r in rows}) == len(rows), 'Duplicate sample records'
paths = [(ROOT/r['narrowPeak']).resolve() for r in rows if r['narrowPeak'] != '-']
assert len(set(paths)) == len(paths), 'Duplicate canonical peak files'
assert set(r['species'] for r in rows) == set(common), 'Unexpected species set'

def count(row):
    path = ROOT/row['narrowPeak']
    q = qc.get(row['sample'], {})
    result = dict(species=row['species'], scientific_name=row['species'].replace('_',' '), common_name=common[row['species']],
                  sample=row['sample'], assembly=row['assembly'], peak_count=None, status='Missing peak file',
                  tss_enrichment=float(q['tss_enrichment']) if q.get('tss_enrichment') not in (None,'','NA','-') else None,
                  final_read_pairs=int(q['final_read_pairs']) if q.get('final_read_pairs') not in (None,'','NA','-') else None,
                  peak_file=row['narrowPeak'], sha256=None, invalid_records=0, duplicate_intervals=0)
    if not path.is_file():
        return result
    n = invalid = duplicates = 0
    seen = set()
    sha = hashlib.sha256()
    with (gzip.open(path, 'rb') if path.suffix == '.gz' else path.open('rb')) as f:
        for lineno, line in enumerate(f, 1):
            sha.update(line)
            if not line.strip() or line.startswith((b'#', b'track', b'browser')):
                continue
            fields = line.split()
            try:
                assert len(fields) >= 3
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                assert start >= 0 and end > start
            except (ValueError, AssertionError):
                invalid += 1
                continue
            key = (chrom, start, end)
            duplicates += key in seen
            seen.add(key)
            n += 1
    assert invalid == 0, f'{path}: {invalid} invalid records'
    result.update(peak_count=n, status='Counted', sha256=sha.hexdigest(), invalid_records=invalid, duplicate_intervals=duplicates)
    return result

if '--reuse-counts' in sys.argv:
    cached = json.loads((OUT/'report_data.json').read_text())
    assert cached['provenance']['manifest_sha256'] == hashlib.sha256(MANIFEST.read_bytes()).hexdigest()
    samples = cached['samples']
else:
    with ThreadPoolExecutor(max_workers=3) as pool:
        samples = list(pool.map(count, rows))
for r in samples:
    r['unique_peak_regions'] = r['peak_count'] - r['duplicate_intervals'] if r['peak_count'] is not None else None
samples.sort(key=lambda r: (r['common_name'], r['sample'].lower()))
groups = defaultdict(list)
for r in samples:
    groups[r['species']].append(r)
summary = []
for species, group in groups.items():
    counts = [r['peak_count'] for r in group if r['peak_count'] is not None]
    unique = [r['unique_peak_regions'] for r in group if r['unique_peak_regions'] is not None]
    summary.append(dict(species=species, scientific_name=species.replace('_',' '), common_name=common[species],
                        samples_listed=len(group), samples_counted=len(counts), missing_peak_files=len(group)-len(counts),
                        median_peaks=statistics.median(counts) if counts else None,
                        mean_peaks=statistics.mean(counts) if counts else None,
                        min_peaks=min(counts) if counts else None, max_peaks=max(counts) if counts else None,
                        sum_sample_peak_calls=sum(counts) if counts else None,
                        median_unique_regions=statistics.median(unique) if unique else None,
                        min_unique_regions=min(unique) if unique else None,
                        max_unique_regions=max(unique) if unique else None))
summary.sort(key=lambda r: -r['median_peaks'])
def save_tsv(name, records):
    with (OUT/name).open('w', newline='') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=list(records[0]))
        w.writeheader(); w.writerows(records)
save_tsv('species_peak_summary.tsv', summary)
save_tsv('sample_peak_counts.tsv', samples)
report = dict(summary=summary, samples=samples, provenance={
    'manifest': str(MANIFEST.relative_to(ROOT)), 'manifest_sha256': hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
    'date': '2026-10-27', 'excluded_species': sorted(EXCLUDE), 'species_count': len(summary),
    'samples_listed':len(samples), 'samples_counted':sum(r['status']=='Counted' for r in samples),
    'missing_samples':[r['sample'] for r in samples if r['status']!='Counted'],
    'total_sample_peak_calls':sum(r['peak_count'] for r in samples if r['peak_count'] is not None),
    'invalid_records':sum(r['invalid_records'] for r in samples), 'duplicate_intervals':sum(r['duplicate_intervals'] for r in samples),
    'method':'Peak entries count valid narrowPeak records (multiple summits within an interval can produce separate records). Distinct peak regions count unique chromosome/start/end tuples within each sample. Distinct regions are not merged across overlapping intervals or samples. Canonical manifest paths used once per sample. No additional peak significance, chromosome, TSS-enrichment or other QC filter. Missing files are unavailable, not zero. Species medians/ranges are across available sample counts; summed calls are not a merged or consensus peak set.'})
(OUT/'report_data.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['provenance'],indent=2))
print('\nCommon name\tSamples\tMedian\tMinimum\tMaximum')
for r in summary:
    print(f"{r['common_name']}\t{r['samples_counted']}\t{r['median_peaks']:,}\t{r['min_peaks']:,}\t{r['max_peaks']:,}")
