"""Indexed queries of precomputed, unweighted gene-body + promoter bigWig sums."""
import csv
import gzip
import json
import sqlite3
import math
import statistics
from functools import lru_cache
from collections import Counter
from contextlib import closing
from pathlib import Path
from region_availability import annotate_region, body_region, MIN_AVAILABLE_FRACTION
from sample_qc import filter_records
import gene_names

DB = Path(__file__).resolve().parent / 'data/gene_activity.sqlite'
DEFAULT_SOURCE = Path(__file__).resolve().parents[2] / 'output/gene_activity_signac_regions'
PEAK_SIGNAL_WINDOWS = {f'peak_signal{k}kb': k * 1000 for k in (20, 50, 100, 200, 500)}


def build(source=DEFAULT_SOURCE, destination=DB):
    source, destination = Path(source), Path(destination)
    manifest = json.loads((source / 'manifest.json').read_text())
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_suffix('.building.sqlite')
    with closing(sqlite3.connect(temp)) as con:
        con.execute('DROP TABLE IF EXISTS activity')
        con.execute('CREATE TABLE activity (species TEXT, symbol TEXT, symbol_key TEXT, locus TEXT, median REAL, unit TEXT, samples TEXT)')
        for species in manifest['species']:
            sid = species['species']
            with gzip.open(source / (sid + '_gene_activity_matrix.tsv.gz'), 'rt') as handle:
                rows = csv.DictReader(handle, delimiter='\t')
                sample_keys = [k for k in rows.fieldnames if k not in {'locus_id','gene_symbol','species_median','valid_samples','unit'}]
                def values():
                    for row in rows:
                        number = lambda s: None if s in ('NA','') else float(s)
                        yield (sid, row['gene_symbol'], row['gene_symbol'].upper(), row['locus_id'], number(row['species_median']), row['unit'], json.dumps({k:number(row[k]) for k in sample_keys}))
                con.executemany('INSERT INTO activity VALUES (?,?,?,?,?,?,?)', values())
        con.execute('CREATE INDEX activity_symbol ON activity(symbol_key)')
        con.execute('DROP TABLE IF EXISTS metadata')
        con.execute('CREATE TABLE metadata (manifest TEXT)')
        con.execute('INSERT INTO metadata VALUES (?)', (json.dumps(manifest),))
        con.commit()
    temp.replace(destination)


def query(gene, strategy="body"):
    if strategy not in ("body", "tss2kb", "tss5kb", "tss100kb", "tss200kb", "tss500kb", "peaks20kb", "peaks50kb", "peaks100kb", "peaks200kb", "peaks500kb") and strategy not in PEAK_SIGNAL_WINDOWS:
        raise ValueError("Unknown gene activity strategy")
    db = DB
    gene = str(gene).strip().upper()
    if not gene or len(gene) > 100:
        raise ValueError('Enter a gene symbol')
    if not db.exists():
        raise ValueError('Gene activity index is missing. Run python activity.py to build it from calculated scores.')
    with closing(sqlite3.connect(f'file:{db}?mode=ro', uri=True)) as con:
        con.row_factory = sqlite3.Row
        symbols = gene_names.candidates(gene)
        placeholders = ','.join('?' for _ in symbols)
        records = [dict(r) for r in con.execute(f'SELECT * FROM activity WHERE symbol_key IN ({placeholders}) ORDER BY species,locus', symbols)]
        if con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='activity_additions'").fetchone():
            records += [dict(r) for r in con.execute(f'SELECT * FROM activity_additions WHERE symbol_key IN ({placeholders}) ORDER BY species,locus', symbols)]
    records, matches = gene_names.match_records(gene, records)
    records = filter_records(refresh_sample_sources(records))
    counts = Counter(r['species'] for r in records)
    for row in records:
        row['ambiguous'] = counts[row['species']] > 1
    if strategy in PEAK_SIGNAL_WINDOWS:
        from peak_activity import score_records
        records = score_records(records, PEAK_SIGNAL_WINDOWS[strategy])
    elif strategy.startswith('peaks'):
        from regulatory_strength import score_records
        records = score_records(records, {'peaks20kb':20000,'peaks50kb':50000,'peaks100kb':100000,'peaks200kb':200000,'peaks500kb':500000}[strategy])
    elif strategy != 'body':
        flank = {'tss2kb':2000, 'tss5kb':5000, 'tss100kb':100000, 'tss200kb':200000, 'tss500kb':500000}[strategy]
        records = tss_records(records, db.stat().st_mtime_ns, flank)
    else:
        from server import SPECS
        folders = {s[0]: Path(__file__).resolve().parent.parent/s[4] for s in SPECS}
        for row in records:
            chrom, start, end = body_region(row['locus'])
            if row.pop('on_demand', False):
                score_record(row, chrom, start, end, db.stat().st_mtime_ns)
            annotate_region(row, chrom, start, end, folders[row['species']])
    lifespans = json.loads((DB.parent / 'species_mls.json').read_text())
    return dict(gene=gene_names.canonical(gene), requested_gene=gene, gene_matches=matches, strategy=strategy, records=records, lifespans=lifespans,
                minimum_available_fraction=MIN_AVAILABLE_FRACTION)


def refresh_sample_sources(records):
    """Use current catalog samples when a stored matrix describes an older set.

    Keep the historical database intact. Replaced/added samples are quantified
    on demand; never reuse or relabel scores from a retired source file.
    """
    from server import SPECS, ROOT, activity_unit
    from sample_qc import EXCLUDED_SAMPLES, sample_id
    present = {r['species'] for r in records}
    sources = {sid: sorted(p for p in (ROOT / folder).glob('*.bw')
                           if p.name.lower().startswith(prefix.lower())
                           and sample_id(p) not in EXCLUDED_SAMPLES)
               for sid, _, _, _, folder, prefix in SPECS if sid in present}
    result = []
    for original in records:
        row = dict(original)
        if row['species'] in sources:
            paths = sources[row['species']]
            samples = json.loads(row['samples']) if isinstance(row['samples'], str) else row['samples']
            retained = {s for s in samples if sample_id(s) not in EXCLUDED_SAMPLES}
            if retained != {p.stem for p in paths}:
                row.update(samples={p.stem: None for p in paths}, median=None, on_demand=True)
                if paths:
                    row['unit'] = activity_unit(paths)
                row.pop('window_details', None)
        result.append(row)
    return result


def score_record(row, chrom, start, end, version):
    from server import SPECS
    folders = {s[0]: Path(__file__).resolve().parent.parent/s[4] for s in SPECS}
    samples, details = {}, {}
    for sample in row['samples']:
        path = folders[row['species']] / (sample+'.bw')
        if not path.exists():
            samples[sample] = None
            details[sample] = {'status':'file_missing','quantified_bp':0}
            continue
        stat = path.stat()
        score, width, status = region_score(str(path),stat.st_size,stat.st_mtime_ns,version,chrom,start,end)
        samples[sample] = score
        details[sample] = {'status':status,'quantified_bp':width}
    valid = [v for v in samples.values() if v is not None]
    row.update(samples=samples, median=statistics.median(valid) if valid else None,
               window_details=details,region_start=start,region_end=end)


@lru_cache(maxsize=8192)
def region_score(path, size, modified, version, chrom, start, end):
    from calculate_gene_activity import quantify
    values, coverage, widths, statuses = quantify(path,[dict(chrom=chrom,region_start=start,region_end=end)])
    value = float(values[0])
    return (value if math.isfinite(value) else None), int(widths[0]), statuses[0]


def register_additional_species():
    """Index loci and sample paths; query bigWigs on demand for new species."""
    from server import SPECS, MODEL_INDEX, EXCLUDED_SAMPLES
    from calculate_gene_activity import gene_regions
    root = Path(__file__).resolve().parent.parent
    with closing(sqlite3.connect(DB)) as con:
        existing = {r[0] for r in con.execute('SELECT DISTINCT species FROM activity')}
        con.execute('DROP TABLE IF EXISTS activity_additions')
        con.execute('CREATE TABLE activity_additions (species TEXT, symbol TEXT, symbol_key TEXT, locus TEXT, median REAL, unit TEXT, samples TEXT, on_demand INTEGER)')
        for sid, name, latin, color, folder, prefix in SPECS:
            if sid in existing:
                continue
            samples = {p.stem:None for p in sorted((root/folder).glob('*.bw'))
                       if p.name.lower().startswith(prefix.lower()) and p.stem.removesuffix('.cpm').removesuffix('_scaled') not in EXCLUDED_SAMPLES}
            if not samples:
                continue
            regions = gene_regions(root/'index'/f'{MODEL_INDEX.get(sid,sid)}_models.tsv.gz')
            sample_json = json.dumps(samples)
            con.executemany('INSERT INTO activity_additions VALUES (?,?,?,?,?,?,?,1)',
                            ((sid,r['gene_symbol'],r['symbol_key'],r['locus_id'],None,'CPM_signal_bp',sample_json) for r in regions))
            print(f'{sid}: {len(regions)} loci indexed for on-demand activity',flush=True)
        con.execute('CREATE INDEX additions_symbol ON activity_additions(symbol_key)')
        con.commit()


def tss_records(records, version, flank=200000):
    from server import ROOT, SPECS
    folders = {s[0]: ROOT/s[4] for s in SPECS}
    result = []
    for original in filter_records(records):
        row = dict(original)
        chrom, span, strand = row['locus'].rsplit('|',1)[1].rsplit(':',2)
        start, end = map(int,span.split('-'))
        tss = start if strand == '+' else end-1
        samples, details = {}, {}
        for sample in row['samples']:
            path = folders[row['species']] / (sample+'.bw')
            if not path.exists():
                samples[sample] = None
                details[sample] = {'status':'file_missing','quantified_bp':0}
                continue
            stat = path.stat()
            score, width, status = window_score(str(path),stat.st_size,stat.st_mtime_ns,version,chrom,tss,flank)
            samples[sample] = score
            details[sample] = {'status':status,'quantified_bp':width}
        valid = [v for v in samples.values() if v is not None]
        row.update(samples=samples, median=statistics.median(valid) if valid else None,
                   region_start=tss-flank,region_end=tss+flank,window_details=details)
        annotate_region(row, chrom, tss-flank, tss+flank, folders[row['species']])
        result.append(row)
    return result


@lru_cache(maxsize=4096)
def window_score(path, size, modified, index_version, chrom, tss, flank=200000):
    from calculate_gene_activity import quantify
    values, coverage, widths, statuses = quantify(path,[dict(chrom=chrom,region_start=tss-flank,region_end=tss+flank)])
    value = float(values[0])
    return (value if math.isfinite(value) else None), int(widths[0]), statuses[0]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination", type=Path, default=DB)
    parser.add_argument('--register-additions', action='store_true', help='Add loci for new catalog species without recomputing existing scores')
    args = parser.parse_args()
    if args.register_additions:
        register_additional_species()
    else:
        build(args.source, args.destination)
    print(args.destination)
