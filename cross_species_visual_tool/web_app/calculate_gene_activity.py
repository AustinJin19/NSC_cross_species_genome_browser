#!/usr/bin/env python3
"""Unweighted bigWig integration over gene bodies + strand-aware promoters.

Uses Signac-style regions, not Signac fragment counts. All annotation biotypes
and gene widths are retained. Never use this matrix as an orthology mapping.
"""
import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict, Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / '.deps'))
import numpy as np
import pyBigWig
from server import SPECS, MODEL_INDEX, EXCLUDED_SAMPLES, resolve_chrom

ROOT = Path(__file__).resolve().parents[2]


def gene_regions(path, upstream=2000, tss_flank=None):
    # The portable browser indexes retain symbols, not universally stable gene IDs.
    # Merge overlapping isoform bodies for a symbol/chromosome/strand; keep
    # disjoint loci separate instead of bridging unrelated duplicated loci.
    groups = defaultdict(list)
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            chrom, start, end, tx, symbol, strand, *_ = line.rstrip().split('\t')
            a, b = int(start), int(end)
            if b > a and strand in ('+', '-'):
                groups[(symbol, chrom, strand)].append((a, b, tx))
    rows = []
    for (symbol, chrom, strand), transcripts in sorted(groups.items()):
        clusters = []
        for a, b, tx in sorted(transcripts):
            if clusters and a < clusters[-1][1]:
                clusters[-1][1] = max(clusters[-1][1], b)
                clusters[-1][2].add(tx)
            else:
                clusters.append([a, b, {tx}])
        for a, b, txs in clusters:
            start, end = (a-upstream, b) if strand == '+' else (a, b+upstream)
            if tss_flank is not None:
                # Outer locus TSS: left edge on +, last body base on -.
                tss = a if strand == '+' else b-1
                start, end = tss-tss_flank, tss+tss_flank
            rows.append(dict(locus_id=f'{symbol}|{chrom}:{a}-{b}:{strand}',
                             gene_symbol=symbol, symbol_key=symbol.upper(), chrom=chrom,
                             strand=strand, body_start=a, body_end=b,
                             region_start=start, region_end=end,
                             transcript_count=len(txs)))
    return sorted(rows, key=lambda r: (r['chrom'], r['body_start'], r['gene_symbol']))


def integrate(intervals, starts, ends):
    """Exact piecewise-constant integrals; unrepresented bases contribute zero.

    Return signal-base integral and number of explicitly represented finite bases.
    Missing chromosomes are handled separately, never converted to zero.
    """
    if not intervals:
        return np.zeros(len(starts)), np.zeros(len(starts))
    data = np.asarray(intervals, dtype=np.float64)
    if not np.all(np.isfinite(data[:, 2])):
        data = data[np.isfinite(data[:, 2])]
    if not len(data):
        return np.zeros(len(starts)), np.zeros(len(starts))
    a, b, v = data.T
    if np.any(a[1:] < b[:-1]):
        raise ValueError('Overlapping bigWig records cannot be integrated as coverage')
    lengths = b-a
    signal = np.r_[0., np.cumsum(lengths*v)]
    coverage = np.r_[0., np.cumsum(lengths)]
    def prefix(x, sums, values):
        idx = np.searchsorted(b, x, side='right')
        total = sums[idx].copy()
        mask = idx < len(a)
        j = idx[mask]
        total[mask] += np.maximum(0, x[mask]-a[j])*values[j]
        return total
    return (prefix(ends, signal, v)-prefix(starts, signal, v),
            prefix(ends, coverage, np.ones(len(a)))-prefix(starts, coverage, np.ones(len(a))))


def quantify(path, rows, chunked=False):
    scores = np.full(len(rows), np.nan)
    represented = np.full(len(rows), np.nan)
    widths = np.zeros(len(rows), dtype=int)
    statuses = ['chromosome_missing']*len(rows)
    with pyBigWig.open(str(path)) as bw:
        sizes = bw.chroms()
        groups = defaultdict(list)
        for i, row in enumerate(rows):
            chrom = resolve_chrom(row['chrom'], sizes)
            if chrom:
                groups[chrom].append(i)
        for chrom, ids in groups.items():
            starts = np.array([max(0, rows[i]['region_start']) for i in ids])
            ends = np.array([min(sizes[chrom], rows[i]['region_end']) for i in ids])
            valid = ends > starts
            keep = np.asarray(ids)[valid]
            for i in np.asarray(ids)[~valid]:
                statuses[i] = 'outside_chromosome'
            qs, qe = starts[valid], ends[valid]
            vals, cov = np.zeros(len(qs)), np.zeros(len(qs))
            if chunked and len(qs):
                # Read each occupied 5 Mb block once; integrate exact stored
                # intervals over every overlapping gene window in that block.
                block = 5_000_000
                for lo in range(int(qs.min())//block*block, int(qe.max()), block):
                    hi = min(lo+block, sizes[chrom])
                    hit = np.flatnonzero((qs < hi) & (qe > lo))
                    if not len(hit):
                        continue
                    records = bw.intervals(chrom, lo, hi) or []
                    clipped = [(max(lo,a),min(hi,b),v) for a,b,v in records]
                    sums, covered = integrate(clipped, np.maximum(qs[hit],lo), np.minimum(qe[hit],hi))
                    vals[hit] += sums
                    cov[hit] += covered
            else:
                for j, (start, end) in enumerate(zip(qs, qe)):
                    total = bw.stats(chrom, int(start), int(end), type='sum', exact=True)[0]
                    fraction = bw.stats(chrom, int(start), int(end), type='coverage', exact=True)[0]
                    vals[j] = total if total is not None else 0.
                    cov[j] = (fraction if fraction is not None else 0.)*(end-start)
            scores[keep], represented[keep], widths[keep] = vals, cov, qe-qs
            for i in keep:
                clipped = rows[i]['region_start'] < 0 or rows[i]['region_end'] > sizes[chrom]
                statuses[i] = 'clipped' if clipped else 'ok'
    return scores, represented, widths, statuses


def value(x):
    return 'NA' if not np.isfinite(x) else f'{x:.10g}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'output/gene_activity_signac_regions')
    parser.add_argument('--upstream', type=int, default=2000)
    parser.add_argument('--tss-flank', type=int, help='Use symmetric fixed TSS windows instead of gene bodies')
    parser.add_argument('--exclude-species', nargs='*', default=[])
    args = parser.parse_args()
    if args.tss_flank is not None and args.tss_flank <= 0:
        parser.error('--tss-flank must be positive')
    if args.upstream < 0:
        parser.error('--upstream must be nonnegative')
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = []
    symbol_scores = defaultdict(dict)
    summary = []
    selected_specs = [s for s in SPECS if s[0] not in args.exclude_species]
    for sid, name, latin, color, folder, prefix in selected_specs:
        annotation = ROOT/'cross_species_visual_tool/index'/f'{MODEL_INDEX.get(sid,sid)}_models.tsv.gz'
        rows = gene_regions(annotation, args.upstream, args.tss_flank)
        samples = sorted(p for p in (ROOT/'cross_species_visual_tool'/folder).glob('*.bw')
                         if p.name.lower().startswith(prefix.lower()) and
                         p.stem.removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES)
        if not samples:
            raise RuntimeError(f'No samples for {sid}')
        scores = []
        unit = 'source_scaled_signal_bp' if sid == 'Human' else 'CPM_signal_bp'
        with gzip.open(args.output/f'{sid}_gene_activity_long.tsv.gz', 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['species','sample',*rows[0].keys(),'signal_sum','signal_mean_per_bp',
                             'quantified_bp','represented_bp','represented_fraction','status','unit'])
            for path in samples:
                v, cov, widths, statuses = quantify(path, rows, chunked=args.tss_flank is not None)
                scores.append(v)
                for i, row in enumerate(rows):
                    w = widths[i]
                    writer.writerow([sid,path.stem,*row.values(),value(v[i]),value(v[i]/w) if w else 'NA',
                                     w,value(cov[i]),value(cov[i]/w) if w else 'NA',statuses[i],unit])
                stat = path.stat()
                provenance.append(dict(species=sid,path=str(path.resolve()),bytes=stat.st_size,
                                       modified_ns=stat.st_mtime_ns,unit=unit,status_counts=dict(Counter(statuses))))
                print(f'{sid}: {path.name}: {len(rows):,} loci quantified', flush=True)
        matrix = np.asarray(scores).T
        medians = np.array([np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else np.nan for v in matrix])
        with gzip.open(args.output/f'{sid}_gene_activity_matrix.tsv.gz','wt') as handle:
            writer=csv.writer(handle,delimiter='\t')
            writer.writerow(['locus_id','gene_symbol',*[p.stem for p in samples],'species_median','valid_samples','unit'])
            for row,v,median in zip(rows,matrix,medians):
                writer.writerow([row['locus_id'],row['gene_symbol'],*[value(x) for x in v],value(median),np.isfinite(v).sum(),unit])
        counts = Counter(row['symbol_key'] for row in rows)
        for row,median in zip(rows,medians):
            symbol_scores[row['symbol_key']][sid] = value(median) if counts[row['symbol_key']]==1 else 'NA'
        summary.append(dict(species=sid,name=name,samples=len(samples),loci=len(rows),symbols=len(counts),
                            ambiguous_symbols=sum(n>1 for n in counts.values()),
                            quantified_loci=int(np.isfinite(matrix).any(axis=1).sum()),unit=unit,
                            annotation=str(annotation),annotation_sha256=hashlib.sha256(annotation.read_bytes()).hexdigest()))
    with (args.output/'species_median_by_symbol.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t');species=[s[0] for s in selected_specs]
        writer.writerow(['symbol_key',*species])
        for symbol,values in sorted(symbol_scores.items()):
            writer.writerow([symbol,*[values.get(s,'NA') for s in species]])
    (args.output/'manifest.json').write_text(json.dumps(dict(region_strategy="tss_window" if args.tss_flank else "gene_body_promoter", tss_flank_bp=args.tss_flank, upstream_bp=args.tss_flank or args.upstream,downstream_bp=args.tss_flank or 0,
        distance_weighting=False,biotypes='all available',max_width=None,
        method='integrated bigWig signal, not fragment counts',species=summary,samples=provenance),indent=2)+'\n')
    with (args.output/'summary.tsv').open('w') as handle:
        keys=['species','name','samples','loci','symbols','ambiguous_symbols','quantified_loci','unit']
        writer=csv.DictWriter(handle,keys,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(summary)


if __name__ == '__main__':
    main()
