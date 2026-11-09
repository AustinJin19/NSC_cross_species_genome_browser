#!/usr/bin/env python3
"""Reproducible five-region activity quantification using current browser inputs.

Run from the project root with PYTHONPATH=cross_species_visual_tool/web_app/.deps
/usr/bin/python3 output/activity_omics_comparison_2026-10-20/quantify_activity.py
"""
import csv
import gzip
import hashlib
import json
import sqlite3
import sys
import time
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'cross_species_visual_tool/web_app').is_dir())
WEB = ROOT / 'cross_species_visual_tool/web_app'
sys.path[:0] = [str(WEB), str(WEB / '.deps')]
import numpy as np
from server import SPECS, EXCLUDED_SAMPLES
from calculate_gene_activity import quantify

OUT = Path(__file__).resolve().parent
CACHE = OUT / 'activity_cache'
CACHE.mkdir(exist_ok=True)
METHODS = [('tss2kb', 2000), ('tss5kb', 5000), ('body', None),
           ('tss200kb', 200000), ('tss500kb', 500000)]
LIFE = json.loads((WEB / 'data/species_mls.json').read_text())
SELECTED = [s for s in SPECS if s[0] not in ('Human', 'Macaque', 'Balaena_mysticetus')
            and LIFE.get(s[2], {}).get('years') and LIFE[s[2]].get('body_mass_g')]


def read_loci(sid):
    with sqlite3.connect(f'file:{WEB}/data/gene_activity.sqlite?mode=ro', uri=True) as con:
        records = []
        for table in ('activity', 'activity_additions'):
            records.extend(con.execute(f'SELECT symbol_key,locus FROM {table} WHERE species=?', (sid,)))
    return sorted(set(records), key=lambda r: r[1])


def regions_for(loci):
    regions = []
    for method, flank in METHODS:
        for gene, locus in loci:
            chrom, span, strand = locus.rsplit('|', 1)[1].rsplit(':', 2)
            a, b = map(int, span.split('-'))
            if flank is None:
                start, end = (a - 2000, b) if strand == '+' else (a, b + 2000)
            else:
                tss = a if strand == '+' else b - 1
                start, end = tss - flank, tss + flank
            regions.append(dict(chrom=chrom, region_start=start, region_end=end,
                                gene=gene, locus=locus, length_kb=(b-a)/1000, method=method))
    return regions


def run_species(spec):
    sid, common, latin, color, folder, prefix = spec
    started = time.time()
    loci = read_loci(sid)
    regions = regions_for(loci)
    files = sorted(p for p in (ROOT / 'cross_species_visual_tool' / folder).glob('*.bw')
                   if p.name.lower().startswith(prefix.lower()) and
                   p.stem.removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES)
    assert files and loci, sid
    fingerprint = dict(version=1, loci_sha256=hashlib.sha256(json.dumps(loci).encode()).hexdigest(),
                       methods=METHODS, files=[dict(path=str(p.resolve()), size=p.stat().st_size,
                                                   mtime_ns=p.stat().st_mtime_ns) for p in files])
    fp = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
    metadata_path = CACHE / f'{sid}.json'
    output_path = CACHE / f'{sid}.tsv.gz'
    if output_path.exists() and metadata_path.exists() and json.loads(metadata_path.read_text()).get('fingerprint') == fp:
        print(f'{sid}: existing verified cache', flush=True)
        return sid
    arrays, fractions = [], []
    requested = np.array([r['region_end']-r['region_start'] for r in regions])
    validation = []
    for path in files:
        sample_cache = CACHE / f'{sid}__{path.stem}.npz'
        if sample_cache.exists():
            cached = np.load(sample_cache)
            if str(cached['fingerprint']) == fp:
                v, w = cached['signal'], cached['widths']
            else:
                v, _, w, _ = quantify(path, regions, chunked=True)
                np.savez_compressed(sample_cache, fingerprint=fp, signal=v, widths=w)
        else:
            v, _, w, _ = quantify(path, regions, chunked=True)
            np.savez_compressed(sample_cache, fingerprint=fp, signal=v, widths=w)
        # Independently check piecewise integration against direct bigWig queries.
        rng = np.random.default_rng(22)
        ids = rng.choice(len(regions), min(35, len(regions)), replace=False)
        direct, _, dw, _ = quantify(path, [regions[i] for i in ids], chunked=False)
        assert np.array_equal(w[ids], dw), (sid, 'interval widths')
        assert np.allclose(v[ids], direct, rtol=2e-7, atol=1e-5, equal_nan=True), (sid, 'integrals', v[ids], direct)
        finite = np.isfinite(v[ids]) & np.isfinite(direct)
        validation.append(dict(sample=path.stem, intervals=len(ids),
                               max_abs_difference=float(np.max(np.abs(v[ids][finite]-direct[finite]))) if finite.any() else 0))
        arrays.append(v)
        fractions.append(np.where(np.isfinite(v), w/requested, np.nan))
        print(f'{sid}: {path.name}; {len(loci):,} loci × 5 regions; {time.time()-started:.1f}s', flush=True)
    matrix = np.asarray(arrays)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        medians = np.nanmedian(matrix, axis=0)
        minimum_fraction = np.nanmin(np.asarray(fractions), axis=0)
    n = np.sum(np.isfinite(matrix), axis=0)
    chosen = {}
    for i, row in enumerate(regions):
        value = medians[i]
        if not np.isfinite(value):
            continue
        key = row['method'], row['gene']
        # Loci already ordered lexically; tie-breaking matches browser queries.
        if key not in chosen or value > medians[chosen[key]]:
            chosen[key] = i
    with gzip.open(output_path, 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['method', 'gene', 'species', 'scientific_name', 'locus', 'median_activity',
                         'length_kb', 'available_fraction', 'region_eligible', 'n_samples'])
        for (method, gene), i in sorted(chosen.items()):
            r = regions[i]
            writer.writerow([method, gene, sid, latin, r['locus'], medians[i], r['length_kb'],
                             minimum_fraction[i], int(minimum_fraction[i] >= .5), n[i]])
    metadata_path.write_text(json.dumps(dict(fingerprint=fp, inputs=fingerprint, species=sid, common=common,
                                             scientific_name=latin, loci=len(loci), validation=validation,
                                             elapsed_seconds=time.time()-started), indent=2))
    print(f'{sid}: COMPLETE in {time.time()-started:.1f}s', flush=True)
    return sid


if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(run_species, spec): spec[0] for spec in SELECTED}
        for future in as_completed(futures):
            print('Finished', future.result(), flush=True)
