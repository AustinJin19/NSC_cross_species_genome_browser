#!/usr/bin/env python3
"""Validate and install species_summary candidate SE shared and union BEDs.

Copies the supplied files unchanged after checking sample IDs, interval bounds
against every contributing sample's bigWig, counts, and exact sample support.
Union (>=1 sample) and shared (>=2 samples) remain separate region sets.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import pyBigWig
import server
from atac_super_enhancers import DATA

ROOT = Path(__file__).resolve().parents[2]


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_bed(path):
    rows = []
    for line in path.read_text().splitlines():
        chrom, start, end, name = line.split('\t')
        start, end = int(start), int(end)
        if start < 0 or end <= start:
            raise ValueError(f'{path.name}: invalid interval {line}')
        rows.append((chrom, start, end, name))
    if len({r[3] for r in rows}) != len(rows):
        raise ValueError(f'{path.name}: duplicate region names')
    return rows


def supported_spans(path, samples, minimum=2):
    """Reconstruct contiguous bases meeting the minimum distinct-sample support."""
    spans = []
    for row in table(path):
        names = row['samples'].split(';')
        if len(set(names)) != int(row['n_samples']) or not set(names).issubset(samples):
            raise ValueError(f'{path.name}: inconsistent sample support')
        if int(row['n_samples']) >= minimum:
            spans.append((row['chrom'], int(row['start']), int(row['end'])))
    merged = []
    for chrom, start, end in sorted(spans):
        if start < 0 or end <= start:
            raise ValueError(f'{path.name}: invalid support interval')
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def install(source, destination):
    server.load_catalog()
    species = table(source / 'species_counts.tsv')
    sample_rows = table(source / 'sample_counts.tsv')
    if len({r['species'] for r in species}) != len(species):
        raise ValueError('Duplicate species in species_counts.tsv')
    manifest = dict(importedAt=datetime.now(timezone.utc).isoformat(),
                    source='species_summary', coordinates='BED: 0-based, half-open',
                    definition={'shared_ge2': 'Bases in candidate super-enhancers from at least two distinct samples; adjacent qualifying segments merged.',
                                'union': 'Union of candidate super-enhancer bases from at least one sample; adjacent qualifying segments merged.'},
                    validation='Shared and union BED/count/support agreement; sample identities and chromosome names/bounds checked against every contributing browser bigWig. No liftOver applied.',
                    sourceTableSha256={name: sha256(source / name) for name in ('species_counts.tsv', 'sample_counts.tsv')},
                    species={})
    for row in species:
        sid = row['species']
        if sid not in server.SAMPLES:
            raise ValueError(f'Unknown browser species: {sid}')
        samples = sorted(r['sample'] for r in sample_rows if r['species'] == sid)
        if len(samples) != int(row['n_samples']) or len(samples) != len(set(samples)):
            raise ValueError(f'{sid}: inconsistent sample count')
        path = source / f'{sid}_candidate_se_shared_ge2.bed'
        support = source / f'{sid}_support_segments.tsv'
        beds = read_bed(path)
        bp = sum(b-a for c,a,b,n in beds)
        if len(beds) != int(row['shared_ge2_regions']) or bp != int(row['shared_ge2_bp']):
            raise ValueError(f'{sid}: BED/count table mismatch')
        if sorted((c,a,b) for c,a,b,n in beds) != supported_spans(support, set(samples)):
            raise ValueError(f'{sid}: BED does not equal the bases with >=2 sample support')
        expected_status = 'single_sample_only' if len(samples) < 2 else 'available'
        if row['replicate_support_status'] != expected_status or (len(samples) < 2 and beds):
            raise ValueError(f'{sid}: inconsistent replicate support status')
        union_path = source / f'{sid}_candidate_se_union.bed'
        union_beds = read_bed(union_path)
        union_bp = sum(b-a for c,a,b,n in union_beds)
        if len(union_beds) != int(row['union_regions']) or union_bp != int(row['union_bp']):
            raise ValueError(f'{sid}: union BED/count table mismatch')
        if sorted((c,a,b) for c,a,b,n in union_beds) != supported_spans(support, set(samples), minimum=1):
            raise ValueError(f'{sid}: union BED does not equal the bases with >=1 sample support')
        for sample in samples:
            if sample not in server.SAMPLES[sid]:
                raise ValueError(f'{sid}: sample {sample} is not in the browser')
            with pyBigWig.open(str(server.SAMPLES[sid][sample])) as bw:
                sizes = bw.chroms()
            if any(c not in sizes or b > sizes.get(c, 0) for c,a,b,n in beds + union_beds):
                raise ValueError(f'{sid}: BED exceeds or mismatches {sample} chromosome bounds')
        manifest['species'][sid] = dict(file=path.name, sha256=sha256(path),
                                        supportSha256=sha256(support), nSamples=len(samples), samples=samples,
                                        regions=len(beds), sharedBp=bp, supportStatus=expected_status,
                                        union=dict(file=union_path.name, sha256=sha256(union_path),
                                                   regions=len(union_beds), coveredBp=union_bp, supportStatus='available'))
    manifest['missingSpecies'] = sorted(set(server.SAMPLES) - set(manifest['species']))
    # Validate everything before installing any file; preserve source BED bytes.
    destination.mkdir(parents=True, exist_ok=True)
    for info in manifest['species'].values():
        shutil.copyfile(source / info['file'], destination / info['file'])
        shutil.copyfile(source / info['union']['file'], destination / info['union']['file'])
    (destination / 'metadata.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(species=len(species), nonemptySpecies=sum(r['regions'] > 0 for r in manifest['species'].values()),
                         sharedRegions=sum(r['regions'] for r in manifest['species'].values()),
                         unionRegions=sum(r['union']['regions'] for r in manifest['species'].values()),
                         missingSpecies=manifest['missingSpecies'], destination=str(destination)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'species_summary')
    parser.add_argument('--destination', type=Path, default=DATA)
    args = parser.parse_args()
    install(args.source, args.destination)
