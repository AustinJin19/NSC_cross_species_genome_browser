#!/usr/bin/env python3
"""Convert GSE167200 Liver top300K contact pixels from hg19 to hg38 BEDPE."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = PROJECT / 'data/interactions/human_liver_deeploop'


def parse_anchor(text):
    """Source bins are 1-based inclusive; return a 0-based half-open BED interval."""
    chrom, interval = text.rsplit(':', 1)
    start, end = map(int, interval.split('-'))
    if start < 1 or end < start:
        raise ValueError(f'Invalid source anchor: {text}')
    return chrom, start - 1, end


def read_source(path):
    rows, anchors = [], {}
    with gzip.open(path, 'rt') as handle:
        for rank, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 4:
                raise ValueError(f'Expected four columns at line {rank}')
            a, b, score, gap = fields
            first, second = parse_anchor(a), parse_anchor(b)
            if first[0] != second[0] or second[1] + 1 - first[2] != int(gap):
                raise ValueError(f'Unexpected source anchor order/gap at line {rank}')
            if not math.isfinite(float(score)):
                raise ValueError(f'Non-finite score at line {rank}')
            for anchor in (a, b):
                if anchor not in anchors:
                    anchors[anchor] = f'a{len(anchors):06d}'
            rows.append((rank, a, b, score, gap))
    return rows, anchors


def valid_mapping(candidates, sizes, source_width):
    if len(candidates) != 1:
        return None, 'unmapped' if not candidates else 'ambiguous'
    chrom, start, end, strand = candidates[0]
    if chrom not in sizes or not 0 <= start < end <= sizes[chrom]:
        return None, 'invalid_target_bounds'
    # Avoid depicting a very large assembly gap as a single precise anchor.
    if not 0.5 <= (end - start) / source_width <= 2:
        return None, 'extreme_length_change'
    return candidates[0], 'mapped'


def sha256(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest() if hasattr(hashlib, 'file_digest') else hashlib.sha256(handle.read()).hexdigest()


def convert(source, chain, liftover, sizes_path, output, min_match=0.95):
    rows, anchors = read_source(source)
    sizes = {f[0]: int(f[1]) for line in sizes_path.read_text().splitlines() if (f := line.split())}
    output.mkdir(parents=True, exist_ok=True)
    candidates = defaultdict(list)
    reasons = {}
    with tempfile.TemporaryDirectory(prefix='human-liver-liftover-') as tmp:
        tmp = Path(tmp)
        source_bed, target_bed, unmapped = (tmp / name for name in ('source.bed', 'hg38.bed', 'unmapped.bed'))
        with source_bed.open('w') as handle:
            for anchor, name in anchors.items():
                chrom, start, end = parse_anchor(anchor)
                handle.write(f'{chrom}\t{start}\t{end}\t{name}\t0\t+\n')
        command = [str(liftover), f'-minMatch={min_match}', str(source_bed), str(chain), str(target_bed), str(unmapped)]
        process = subprocess.run(command, check=True, capture_output=True, text=True)
        for line in target_bed.read_text().splitlines():
            f = line.split('\t')
            candidates[f[3]].append((f[0], int(f[1]), int(f[2]), f[5]))
        reason = 'unmapped'
        for line in unmapped.read_text().splitlines():
            if line.startswith('#'):
                reason = line[1:].strip()
            elif line:
                reasons[line.split('\t')[3]] = reason
    mapping, status = {}, {}
    with gzip.open(output / 'anchor_mapping.tsv.gz', 'wt') as audit:
        audit.write('source_hg19_1based\tchrom_hg38\tstart0_hg38\tend_hg38\tmapping_strand\tstatus\tdetail\n')
        for anchor, name in anchors.items():
            _, start, end = parse_anchor(anchor)
            mapped, flag = valid_mapping(candidates[name], sizes, end - start)
            status[anchor] = flag
            if mapped:
                mapping[anchor] = mapped
            audit.write('\t'.join(map(str, (anchor, *(mapped or ('.', '.', '.', '.')), flag, reasons.get(name, '')))) + '\n')
    bedpe = output / 'GSE167200_Liver.top300K.hg38.bedpe.gz'
    rejected_reasons = Counter()
    retained = 0
    with gzip.open(bedpe, 'wt') as out, gzip.open(output / 'rejected_pairs.tsv.gz', 'wt') as rejected:
        out.write('#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\thg19_anchor1\thg19_anchor2\thg19_gap_bp\tmapping_strand1\tmapping_strand2\n')
        rejected.write('source_rank\thg19_anchor1\thg19_anchor2\tscore\treason\n')
        for rank, a, b, score, gap in rows:
            reason = ''
            if a not in mapping or b not in mapping:
                reason = f'anchor1:{status[a]};anchor2:{status[b]}'
            elif mapping[a][0] != mapping[b][0]:
                reason = 'cis_pair_mapped_to_different_chromosomes'
            if reason:
                rejected_reasons[reason] += 1
                rejected.write(f'{rank}\t{a}\t{b}\t{score}\t{reason}\n')
                continue
            first, second = mapping[a], mapping[b]
            fields = (*first[:3], *second[:3], f'GSE167200_Liver_{rank:06d}', score, '.', '.', a, b, gap, first[3], second[3])
            out.write('\t'.join(map(str, fields)) + '\n')
            retained += 1
    summary = dict(source_file=source.name, source_sha256=sha256(source), source_assembly='hg19', target_assembly='hg38',
                   source_coordinates='1-based inclusive, inferred from contiguous source bins (next start = previous end + 1) and minimum start 1',
                   output_coordinates='BEDPE: 0-based starts, exclusive ends; strand columns are . because contacts are unstranded',
                   chain_file=chain.name, chain_sha256=sha256(chain), min_match=min_match,
                   anchor_length_ratio_limits=[0.5, 2], unique_source_anchors=len(anchors), anchor_status=dict(Counter(status.values())),
                   source_pairs=len(rows), retained_pairs=retained, rejected_pairs=len(rows)-retained, rejected_reasons=dict(rejected_reasons),
                   bedpe_file=bedpe.name, bedpe_sha256=sha256(bedpe),
                   description='Human liver DeepLoop top300K ranked contact pixels, not FDR-filtered independent or experimentally validated loops. Scores are original model-derived loop strengths, not p-values.',
                   source_urls=['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167200',
                                'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167200/suppl/GSE167200_Liver.top300K.txt.gz',
                                'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167200/suppl/GSE167200_README.txt',
                                'https://github.com/JinLabBioinfo/DeepLoop'], liftover_log=process.stderr.strip())
    (output / 'conversion_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=PROJECT / 'GSE167200_Liver.top300K.txt.gz')
    parser.add_argument('--chain', type=Path, default=PROJECT / 'data/chipseq/ChIP_seq/liftover/chains/hg19ToHg38.over.chain.gz')
    parser.add_argument('--liftover', type=Path, default=PROJECT / 'data/chipseq/ChIP_seq/liftover/tools/liftOver')
    parser.add_argument('--sizes', type=Path, default=PROJECT / 'reference_genomes/selected3/hg38.fa.fai')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--min-match', type=float, default=0.95)
    args = parser.parse_args()
    print(json.dumps(convert(args.source, args.chain, args.liftover, args.sizes, args.output, args.min_match), indent=2))
