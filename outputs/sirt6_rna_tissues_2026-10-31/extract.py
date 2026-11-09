"""Read the complete RNA workbooks, retaining locus genes and tissue metadata."""
import argparse
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT/'cross_species_visual_tool/web_app'))
from build_omics import rows, number

parser = argparse.ArgumentParser()
parser.add_argument('workbook')
args = parser.parse_args()
path = ROOT/'data/omics/rna'/args.workbook
genes = json.loads((ROOT/'outputs/sirt6_neighborhood_rank_matrices_2026-10-29/matrix_data.json').read_text())['metadata']['genes']
lookup = {g.upper():g for g in genes}
iterator = rows(path)
_, header = next(iterator)
columns = {v:k for k,v in header.items()}
records = []
samples = {}
counts = Counter()
total = 0
started = time.time()
for row_number, row in iterator:
    total += 1
    get = lambda key: row.get(columns.get(key))
    sample, species = get('sample'), get('species')
    original_tissue = get('tissue')
    # Previously confirmed by the user; the source workbook mislabeled this lung sample.
    tissue = 'Lung' if sample == 'BMR_Turk_Lu' else str(original_tissue or '').strip().title()
    key = (sample, species, tissue)
    if key not in samples:
        samples[key] = dict(sample=sample, species=species, tissue=tissue,
            source_tissue=original_tissue, common=get('common'), source=get('source'), rows=0)
    samples[key]['rows'] += 1
    counts[tissue] += 1
    symbol = str(get('symbol') or '').upper()
    if symbol in lookup:
        records.append(dict(gene=lookup[symbol], source_symbol=get('symbol'), sample=sample,
            species=species, common=get('common'), tissue=tissue, source_tissue=original_tissue,
            tpm=number(get('tpm')), log_tpm=number(get('log.tpm')), counts=number(get('counts')),
            source=get('source'), workbook=path.name, source_row=row_number))
    if total % 250000 == 0:
        print(f'{path.name}: {total:,} rows; {len(records):,} locus records', flush=True)
payload = dict(workbook=path.name, path=str(path), bytes=path.stat().st_size,
    sha256=hashlib.sha256(path.read_bytes()).hexdigest(), total_rows=total,
    tissue_rows=dict(counts), samples=list(samples.values()), records=records)
(OUT/(path.stem+'_extracted.json')).write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
print(json.dumps(dict(file=path.name, rows=total, locus_rows=len(records), tissues=dict(counts),
    samples=len(samples), elapsed_seconds=round(time.time()-started,1))), flush=True)
