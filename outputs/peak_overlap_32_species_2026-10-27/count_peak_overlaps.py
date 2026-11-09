#!/usr/bin/env python3
"""Coordinate overlap within samples and shared coverage within each species.

BED coordinates are zero-based, end-exclusive. Cross-sample shared regions are
maximal contiguous segments with coverage by >=2, or all, available samples.
Adjacent qualifying segments are joined, even if the supporting pair changes.
This is coordinate overlap (>=1 bp), not IDR or a reciprocal-overlap filter.
"""
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import groupby, combinations
import csv, gzip, hashlib, heapq, json, random, statistics

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SOURCE = ROOT/'outputs/peak_counts_32_species_2026-10-27/report_data.json'

def merge_intervals(intervals):
    merged = []
    for a,b in sorted(intervals):
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(b,merged[-1][1]))
        else:
            merged.append((a,b))
    return merged

def within_overlap(intervals):
    """Count distinct intervals in positive-overlap components, ignoring duplicates."""
    clusters = members = 0
    end = -1
    size = 0
    for a,b in sorted(intervals):
        if a < end:
            size += 1
            end = max(end,b)
        else:
            if size > 1:
                clusters += 1
                members += size
            size,end = 1,b
    if size > 1:
        clusters += 1
        members += size
    return clusters,members

def shared_coverage(sample_intervals):
    """Return maximal connected segment counts and bp for each support rule."""
    n = len(sample_intervals)
    rules = {'union':lambda m: m != 0}
    if n >= 2:
        rules['at_least_two'] = lambda m: m.bit_count() >= 2
        rules['all_samples'] = lambda m: m == (1<<n)-1
        for i,j in combinations(range(n),2):
            bitmask = (1<<i)|(1<<j)
            rules[f'pair_{i}_{j}'] = lambda m, b=bitmask: m & b == b
    counts = {k:0 for k in rules}
    bp = {k:0 for k in rules}
    def events(intervals, i):
        for a,b in intervals:
            yield a,1,i
            yield b,-1,i
    chroms = sorted(set().union(*(set(d) for d in sample_intervals)))
    for chrom in chroms:
        streams = [events(d.get(chrom,[]),i) for i,d in enumerate(sample_intervals)]
        merged_events = heapq.merge(*streams)
        active = [0]*n
        mask = 0
        last = None
        last_end = {k:None for k in rules}
        for pos,grp in groupby(merged_events, key=lambda x:x[0]):
            if last is not None and pos > last:
                for key,test in rules.items():
                    if test(mask):
                        if last_end[key] != last:
                            counts[key] += 1
                        bp[key] += pos-last
                        last_end[key] = pos
            for _,delta,i in grp:
                active[i] += delta
                assert active[i] in (0,1)
                if active[i]: mask |= 1<<i
                else: mask &= ~(1<<i)
            last = pos
        assert mask == 0
    return {k:{'regions':counts[k],'bp':bp[k]} for k in rules}

def load_sample(row):
    if row['peak_count'] is None:
        return dict(row,duplicate_coordinate_groups=None,distinct_intervals_with_overlap=None,
                    within_sample_overlap_clusters=None,covered_bp=None),None
    path = ROOT/row['peak_file']
    by_chrom = defaultdict(Counter)
    sha = hashlib.sha256()
    raw = 0
    with (gzip.open(path,'rb') if path.suffix=='.gz' else path.open('rb')) as handle:
        for line in handle:
            sha.update(line)
            if not line.strip() or line.startswith((b'#',b'track',b'browser')): continue
            f = line.split()
            a,b = int(f[1]),int(f[2])
            assert 0 <= a < b
            by_chrom[f[0].decode()][(a,b)] += 1
            raw += 1
    unique = sum(len(v) for v in by_chrom.values())
    assert sha.hexdigest() == row['sha256'], f'File changed: {path}'
    assert raw == row['peak_count'] and unique == row['unique_peak_regions']
    clusters = members = duplicate_groups = 0
    merged = {}
    for chrom,counts in by_chrom.items():
        c,m = within_overlap(counts)
        clusters += c; members += m
        duplicate_groups += sum(v>1 for v in counts.values())
        merged[chrom] = merge_intervals(counts)
    result = dict(row,duplicate_coordinate_groups=duplicate_groups,
                  distinct_intervals_with_overlap=members,within_sample_overlap_clusters=clusters,
                  covered_bp=sum(b-a for intervals in merged.values() for a,b in intervals))
    return result,merged

def process_species(item):
    species,rows = item
    rows = sorted(rows,key=lambda r:r['sample'])
    samples,coverage,available = [],[],[]
    for r in rows:
        s,c = load_sample(r)
        samples.append(s)
        if c is not None: coverage.append(c); available.append(s)
    n = len(available)
    assert n > 0
    stats = shared_coverage(coverage)
    assert stats['union']['bp'] <= sum(s['covered_bp'] for s in available)
    if n>=2:
        assert stats['all_samples']['bp'] <= stats['at_least_two']['bp'] <= stats['union']['bp']
    summary = dict(species=species,common_name=rows[0]['common_name'],samples_available=n,
                   samples_listed=len(rows),union_regions=stats['union']['regions'],
                   shared_at_least_two=stats.get('at_least_two',{}).get('regions'),
                   shared_all_samples=stats.get('all_samples',{}).get('regions'),
                   union_bp=stats['union']['bp'],
                   shared_at_least_two_bp=stats.get('at_least_two',{}).get('bp'),
                   shared_all_samples_bp=stats.get('all_samples',{}).get('bp'),
                   median_peak_entries=statistics.median(s['peak_count'] for s in available),
                   median_distinct_intervals=statistics.median(s['unique_peak_regions'] for s in available),
                   median_within_sample_overlapping_distinct_intervals=statistics.median(s['distinct_intervals_with_overlap'] for s in available),
                   median_extra_same_coordinate_entries=statistics.median(s['duplicate_intervals'] for s in available),
                   sample_names='; '.join(s['sample'] for s in available),
                   missing_samples='; '.join(s['sample'] for s in samples if s['peak_count'] is None))
    pairwise=[]
    for i,j in combinations(range(n),2):
        a,b=available[i],available[j]
        inter=stats[f'pair_{i}_{j}']
        union_bp=a['covered_bp']+b['covered_bp']-inter['bp']
        pairwise.append(dict(species=species,common_name=rows[0]['common_name'],sample_a=a['sample'],sample_b=b['sample'],
                             shared_regions=inter['regions'],shared_bp=inter['bp'],
                             covered_bp_a=a['covered_bp'],covered_bp_b=b['covered_bp'],
                             jaccard_bp=inter['bp']/union_bp if union_bp else None))
    return summary,samples,pairwise

def validate():
    # Independent base-resolution check on synthetic data: ties, nested regions,
    # duplicates, book-ended intervals, disjoint contigs, and changing support.
    rng=random.Random(42)
    for n in (1,2,3):
        for _ in range(60):
            maps=[]; bitsets=[]
            for i in range(n):
                raw=[(a,a+rng.randint(1,12)) for a in (rng.randrange(40) for k in range(12))]
                maps.append({'chr1':merge_intervals(raw)})
                bitsets.append(set(x for a,b in raw for x in range(a,b)))
                uniq=sorted(set(raw))
                expected=sum(any(max(a,c)<min(b,d) for c,d in uniq if (a,b)!=(c,d)) for a,b in uniq)
                assert within_overlap(uniq)[1]==expected
            stats=shared_coverage(maps)
            masks={x:sum(1<<i for i,s in enumerate(bitsets) if x in s) for x in range(52)}
            for key,result in stats.items():
                if key=='union': selected={x for x,m in masks.items() if m}
                elif key=='at_least_two': selected={x for x,m in masks.items() if m.bit_count()>=2}
                elif key=='all_samples': selected=set.intersection(*bitsets)
                else:
                    _,i,j=key.split('_'); selected=bitsets[int(i)]&bitsets[int(j)]
                assert result['bp']==len(selected)
                assert result['regions']==sum(x-1 not in selected for x in selected)
    assert shared_coverage([{'a':[(0,5)]},{'a':[(5,10)]}])['all_samples']['regions']==0
    assert shared_coverage([{'a':[(0,10)]},{'b':[(0,10)]}])['all_samples']['regions']==0
    print('Validated sweep against 180 independent base-resolution examples.',flush=True)

def save_tsv(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)

if __name__=='__main__':
    validate()
    prior=json.loads(SOURCE.read_text())
    groups=defaultdict(list)
    for r in prior['samples']: groups[r['species']].append(r)
    summaries=[];samples=[];pairs=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(process_species,item) for item in groups.items()]
        for future in as_completed(futures):
            summary,s,p=future.result()
            summaries.append(summary);samples.extend(s);pairs.extend(p)
            print(f"{len(summaries):2d}/32 {summary['common_name']}: n={summary['samples_available']}, shared >=2={summary['shared_at_least_two']}, all={summary['shared_all_samples']}",flush=True)
    summaries.sort(key=lambda r:(r['shared_at_least_two'] is None,-(r['shared_at_least_two'] or 0),r['common_name']))
    samples.sort(key=lambda r:(r['common_name'],r['sample']))
    pairs.sort(key=lambda r:(r['common_name'],r['sample_a'],r['sample_b']))
    save_tsv('species_peak_overlaps.tsv',summaries)
    save_tsv('sample_peak_overlaps.tsv',samples)
    save_tsv('pairwise_peak_overlaps.tsv',pairs)
    report=dict(species=summaries,samples=samples,pairwise=pairs,method={
        'date':'2026-10-27','excluded_species':prior['provenance']['excluded_species'],
        'input_manifest':prior['provenance']['manifest'],'manifest_sha256':prior['provenance']['manifest_sha256'],
        'sample_hashes_verified_against_prior_counts':True,
        'coordinate_system':'BED zero-based, half-open; chromosome names must match exactly; strand ignored.',
        'shared_region':'Maximal contiguous base-pair segment covered by peaks in at least two, or all, available samples of one species. Require at least 1 bp, no reciprocal-overlap cutoff. Adjacent qualifying segments are joined. Counts can split a broad peak into multiple shared segments.',
        'duplicates':'Deduplicate identical chromosome/start/end intervals within each sample before overlap; per-sample overlapping/book-ended intervals merged for coverage. A sample contributes at most once at each base.',
        'within_sample':'Distinct intervals overlapping another distinct interval by >=1 bp. Exact-coordinate extra entries are reported separately.',
        'single_sample':'Cross-sample overlap unavailable, not zero.',
        'filters':'Same 80 available samples across 32 species. No additional TSS, peak significance, chromosome or mappability filter. Missing NMR1Liver omitted.',
        'interpretation':'Coordinate overlap only; not an IDR reproducibility analysis. Counts depend on peak widths, sample quality, and sample number.'})
    (OUT/'overlap_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Saved species, sample, and pairwise overlap tables.',flush=True)
