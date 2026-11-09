"""Build a portable, evidence-labelled gene-name resolver from HGNC/HCOP.

Run with --sources containing the three human_*_hcop_fifteen_column.txt.gz
files and --hgnc pointing to the HGNC complete-set TSV. No network at runtime.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
HCOP = 'https://www.genenames.org/tools/hcop/'
HGNC = 'https://www.genenames.org/download/statistics-and-files/'


def build(sources, hgnc, destination):
    from server import SPECS, MODEL_INDEX, ROOT
    approved, alias_targets = set(), defaultdict(set)
    with hgnc.open() as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            if r['status'] != 'Approved':
                continue
            canonical = r['symbol'].upper()
            approved.add(canonical)
            for field in ('alias_symbol', 'prev_symbol'):
                for alias in r[field].split('|'):
                    if alias:
                        alias_targets[alias.upper()].add(canonical)
    # Never reinterpret an approved human symbol as another gene's alias.
    aliases = {a: next(iter(v)) for a, v in alias_targets.items()
               if len(v) == 1 and a not in approved}
    orthologs, provenance = {}, []
    for sid, species in [('Mouse', 'mouse'), ('Rat', 'rat'), ('Macaque', 'macaque')]:
        file = sources/f'human_{species}_hcop_fifteen_column.txt.gz'
        pairs = {}
        with gzip.open(file, 'rt') as handle:
            for r in csv.DictReader(handle, delimiter='\t'):
                human, target = r['human_symbol'].upper(), r[f'{species}_symbol']
                support = set(r['support'].split(','))
                if human == '-' or target == '-':
                    continue
                key = (human, target.upper())
                pair = pairs.setdefault(key, dict(symbol=target, support=set(), ids=set()))
                pair['support'].update(support)
                pair['ids'].add(r[f'{species}_entrez_gene'])
        # Require agreement plus a gene-level source; broad family-only calls
        # must not merge paralogs (especially the globin gene clusters).
        pairs = {k: v for k, v in pairs.items() if len(v['support']) >= 3
                 and v['support'] & {'Ensembl', 'NCBI', 'HGNC'}}
        forward, reverse = defaultdict(set), defaultdict(set)
        for human, target in pairs:
            forward[human].add(target)
            reverse[target].add(human)
        mapping = {}
        for human, targets in forward.items():
            entries = []
            for target in sorted(targets):
                p = pairs[(human, target)]
                entries.append(dict(symbol=p['symbol'], support=sorted(p['support']),
                                    target_gene_ids=sorted(p['ids']-{'-'}),
                                    reciprocal_unique=len(targets) == 1 and len(reverse[target]) == 1))
            mapping[human] = entries
        orthologs[sid] = mapping
        provenance.append(dict(file=file.name, sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
                               url=f'https://storage.googleapis.com/public-download-files/hcop/{file.name}'))
    projected = {}
    for sid, _, _, _, _, _ in SPECS:
        names, mouse_names = {}, set()
        with gzip.open(ROOT/'index'/f'{MODEL_INDEX.get(sid, sid)}_models.tsv.gz', 'rt') as handle:
            for line in handle:
                f = line.rstrip().split('\t')
                names.setdefault(f[4].upper(), f[4])
                if f[3].startswith('ENSMUST') and '#' in f[3]:
                    mouse_names.add(f[4].upper())
        if sid != 'Mouse' and mouse_names:
            projected[sid] = sorted(mouse_names)
    payload = dict(version=1, built=str(date.today()), approved=sorted(approved), aliases=aliases,
                   orthologs=orthologs, projected_mouse_symbols=projected,
                   species_latin={s[0]:s[2] for s in SPECS}, sources=provenance+[
                       dict(file=hgnc.name, sha256=hashlib.sha256(hgnc.read_bytes()).hexdigest(), url=HGNC)],
                   policy='HGNC unambiguous aliases; HCOP >=3 sources including Ensembl/NCBI/HGNC; reciprocal-unique pairs only for automatic ortholog matching. ENSMUST-derived annotation names are labelled projected, not independently verified orthology.')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(destination, 'wt') as handle:
        json.dump(payload, handle, separators=(',', ':'))
    print(json.dumps(dict(aliases=len(aliases), ortholog_genes={s:len(g) for s,g in orthologs.items()},
                          projected_annotations={s:len(g) for s,g in projected.items()}), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--hgnc', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE/'data/gene_names.json.gz')
    args = parser.parse_args()
    build(args.sources, args.hgnc, args.output)
