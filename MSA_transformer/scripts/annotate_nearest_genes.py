"""Annotate BED-style regions with GENCODE transcript TSS and gene proximity.

Standard library only. Retains ties, input row order, and original cell strings.
"""
import argparse
import bisect
import csv
import gzip
import hashlib
import json
import random
import re
import statistics
from collections import defaultdict
from pathlib import Path


def distance(start, end, point):
    """Distance between nucleotide positions; region is [start, end)."""
    return max(start - point, point - (end - 1), 0)


def interval_distance(start, end, gs, ge):
    return max(gs - (end - 1), start - (ge - 1), 0)


def nearest(start, end, positions):
    i, j = bisect.bisect_left(positions, start), bisect.bisect_left(positions, end)
    if i < j:
        return 0, positions[i:j]
    candidates = positions[max(0, i - 1): min(len(positions), i + 1)]
    if not candidates:
        return None, []
    best = min(distance(start, end, p) for p in candidates)
    return best, [p for p in candidates if distance(start, end, p) == best]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regions', required=True, type=Path)
    parser.add_argument('--gtf', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    with args.regions.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        original_fields = reader.fieldnames
        rows = list(reader)
    chroms = {r['chrom'] for r in rows}
    genes, tss = {}, defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    attributes = re.compile(r'(\w+) "([^"]*)"')
    with gzip.open(args.gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[0] not in chroms or fields[2] not in ('gene', 'transcript'):
                continue
            chrom, _, feature, start, end, _, strand, _, text = fields
            a = dict(attributes.findall(text))
            start, end = int(start) - 1, int(end)
            gid = a['gene_id']
            if feature == 'gene':
                assert gid not in genes, gid
                genes[gid] = dict(chrom=chrom, start=start, end=end, strand=strand,
                                  symbol=a.get('gene_name', gid), biotype=a['gene_type'])
            else:
                pos = start if strand == '+' else end - 1
                tss[chrom][pos][gid].add(a['transcript_id'])
    by_chrom = defaultdict(list)
    for gid, gene in genes.items():
        by_chrom[gene['chrom']].append(gid)
    positions = {c: sorted(tss[c]) for c in chroms}
    coding_positions = {c: [p for p in positions[c]
        if any(genes[g]['biotype'] == 'protein_coding' for g in tss[c][p])]
        for c in chroms}
    assert all(positions.values()) and all(coding_positions.values())
    annotated, details = [], []
    detail_fields = ['region_row', 'chrom', 'start', 'end', 'candidate_set',
                     'gene_id', 'gene_symbol', 'gene_biotype', 'gene_strand',
                     'tss_0based', 'distance_to_tss_bp', 'region_relative_to_tss',
                     'transcript_ids', 'overlaps_gene_body']

    def sort_genes(gids):
        return sorted(set(gids), key=lambda g: (genes[g]['symbol'], g))

    def add_gene_fields(out, prefix, gids):
        gids = sort_genes(gids)
        out[prefix + '_gene_symbols'] = ';'.join(genes[g]['symbol'] for g in gids) or 'NA'
        out[prefix + '_gene_ids'] = ';'.join(gids) or 'NA'
        out[prefix + '_gene_biotypes'] = ';'.join(genes[g]['biotype'] for g in gids) or 'NA'
        out[prefix + '_gene_count'] = len(gids)

    sample = set(random.Random(42).sample(range(len(rows)), min(100, len(rows))))
    for idx, row in enumerate(rows):
        chrom, start, end = row['chrom'], int(row['start']), int(row['end'])
        assert 0 <= start < end and int(row['Length']) == end - start
        out = dict(row)
        out['region_row'] = idx + 1
        for prefix, pos_list in [('nearest_tss', positions[chrom]),
                                 ('nearest_protein_coding_tss', coding_positions[chrom])]:
            best, hits = nearest(start, end, pos_list)
            # Independent exhaustive check of the indexed search, including ties.
            if idx in sample:
                brute_best = min(distance(start, end, p) for p in pos_list)
                brute_hits = [p for p in pos_list if distance(start, end, p) == brute_best]
                assert (best, hits) == (brute_best, brute_hits)
            gids = set()
            for pos in hits:
                for gid, tids in sorted(tss[chrom][pos].items()):
                    gene = genes[gid]
                    if prefix == 'nearest_protein_coding_tss' and gene['biotype'] != 'protein_coding':
                        continue
                    gids.add(gid)
                    overlaps = start < gene['end'] and end > gene['start']
                    if start <= pos < end:
                        relative = 'overlaps_TSS'
                    else:
                        left = end <= pos
                        relative = 'upstream' if left == (gene['strand'] == '+') else 'downstream'
                    details.append(dict(region_row=idx + 1, chrom=chrom, start=start, end=end,
                        candidate_set=prefix, gene_id=gid, gene_symbol=gene['symbol'],
                        gene_biotype=gene['biotype'], gene_strand=gene['strand'], tss_0based=pos,
                        distance_to_tss_bp=best, region_relative_to_tss=relative,
                        transcript_ids=';'.join(sorted(tids)), overlaps_gene_body=int(overlaps)))
            add_gene_fields(out, prefix, gids)
            out[prefix + '_distance_bp'] = best if best is not None else 'NA'
        body_distances = {g: interval_distance(start, end, genes[g]['start'], genes[g]['end'])
                          for g in by_chrom[chrom]}
        best_body = min(body_distances.values())
        add_gene_fields(out, 'nearest_gene_body', [g for g, d in body_distances.items() if d == best_body])
        out['nearest_gene_body_distance_bp'] = best_body
        add_gene_fields(out, 'overlapping', [g for g, d in body_distances.items() if d == 0])
        annotated.append(out)

    # Exercise one-base intervals, ties, overlap, and edge-of-chromosome lookup.
    assert nearest(100, 101, [99, 101]) == (1, [99, 101])
    assert nearest(100, 105, [99, 100, 104, 105]) == (0, [100, 104])
    assert nearest(0, 1, [1]) == (1, [1])
    assert nearest(110, 111, [100]) == (10, [100])
    assert interval_distance(100, 101, 101, 110) == 1
    assert interval_distance(100, 105, 104, 110) == 0
    args.outdir.mkdir(parents=True, exist_ok=True)
    main_path = args.outdir / (args.regions.stem + '_nearest_genes_hg38.tsv')
    detail_path = args.outdir / 'nearest_TSS_candidates_hg38.tsv'
    for path, data, fields in [(main_path, annotated, list(annotated[0])),
                               (detail_path, details, detail_fields)]:
        with path.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
    with main_path.open() as f:
        saved = list(csv.DictReader(f, delimiter='\t'))
    assert len(saved) == len(rows)
    for src, dest in zip(rows, saved):
        assert all(src[k] == dest[k] for k in original_fields)
        assert dest['nearest_tss_gene_ids'] != 'NA'
        assert dest['nearest_protein_coding_tss_gene_ids'] != 'NA'
    summary = dict(region_count=len(rows), annotation='GENCODE v50 comprehensive CHR, GRCh38.p14',
        input_sha256=sha(args.regions), gtf_sha256=sha(args.gtf),
        coordinate_assumption='Input BED-style 0-based half-open; Length equals end-start.',
        reference_url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz',
        annotated_genes_on_input_chromosomes=len(genes),
        regions_overlapping_gene_bodies=sum(r['overlapping_gene_count'] > 0 for r in annotated),
        nearest_tss_tied_gene_regions=sum(r['nearest_tss_gene_count'] > 1 for r in annotated),
        nearest_coding_tss_tied_gene_regions=sum(r['nearest_protein_coding_tss_gene_count'] > 1 for r in annotated),
        all_original_cells_preserved=True, exhaustive_distance_checks=2 * len(sample),
        tss_candidate_rows=len(details))
    for prefix in ['nearest_tss', 'nearest_protein_coding_tss']:
        ds = [r[prefix + '_distance_bp'] for r in annotated]
        summary[prefix] = dict(median_distance_bp=statistics.median(ds),
                              overlaps_tss=sum(d == 0 for d in ds),
                              within_2000bp=sum(d <= 2000 for d in ds),
                              within_10000bp=sum(d <= 10000 for d in ds),
                              beyond_100000bp=sum(d > 100000 for d in ds))
    (args.outdir / 'nearest_genes_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    print('First five nearest protein-coding TSS assignments:')
    for r in annotated[:5]:
        print(r['chrom'], r['start'], r['end'], r['nearest_protein_coding_tss_gene_symbols'],
              r['nearest_protein_coding_tss_distance_bp'])


if __name__ == '__main__':
    main()
