"""Prepare the original 71 genes from the validated current 500-kb cache."""
import csv
import hashlib
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SUPPORT = Path(__file__).resolve().parent
OUT = SUPPORT.parent
ROOT = OUT.parents[1]
SOURCE = ROOT / 'outputs/regulatory_mls_plot_2026-11-07'
ORIGINAL = ROOT / 'outputs/regulatory-500kb-mass-2026-10-19'
APP = ROOT / 'cross_species_visual_tool/web_app'
sys.path[:0] = [str(APP), str(APP / '.deps')]
import gene_names
import server


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((SOURCE / 'fdr_manifest.json').read_text())
    # Refuse to silently combine changed inputs with the saved BH families.
    for entry in manifest['source_files'] + manifest['metadata_files']:
        path = Path(entry['path'])
        stat = path.stat() if path.exists() else None
        if (stat.st_size if stat else None) != entry['size'] or \
           (stat.st_mtime_ns if stat else None) != entry['mtime_ns']:
            raise ValueError('Rerun the genome-wide scan before plotting: ' + str(path))
    genes = (ORIGINAL / 'FDR05_positive_genes.txt').read_text().split()
    assert len(genes) == len(set(genes)) == 71
    resolver = gene_names.resolver()
    assert all(resolver.canonical(g) == g for g in genes), 'Review changed canonical gene names'
    life = json.loads((APP / 'data/species_mls.json').read_text())
    rows = []
    cache_fingerprints = {}
    for sid, common, latin, *_ in server.SPECS:
        trait = life.get(latin, {})
        if sid in manifest['excluded_species'] or not trait.get('years') or trait.get('lower_bound_years'):
            continue
        provenance = json.loads((SOURCE / 'fdr_quantification' / (sid + '.json')).read_text())
        digest = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
        with np.load(SOURCE / 'fdr_quantification' / (sid + '.npz')) as z:
            assert str(z['fingerprint']) == digest, sid + ': cache does not match provenance'
            medians, available = z['medians'][0], z['available'][0]
        loci = provenance['loci']
        assert len(loci) == len(medians) == len(available)
        by_symbol = defaultdict(list)
        for i, (symbol, locus) in enumerate(loci):
            by_symbol[symbol].append(i)
        cache_fingerprints[sid] = digest
        for gene in genes:
            match = resolver.resolve(gene, sid, by_symbol)
            indexes = [i for i in by_symbol.get(match['symbol'], []) if np.isfinite(medians[i])]
            i = min(indexes, key=lambda j: (-medians[j], loci[j][1])) if indexes else None
            median = float(medians[i]) if i is not None else None
            fraction = float(available[i]) if i is not None and np.isfinite(available[i]) else None
            rows.append(dict(gene=gene, strategy='peaks500kb', species=sid,
                common_name=common, scientific_name=latin,
                matched_symbol=match['symbol'], locus=loci[i][1] if i is not None else '',
                activity=median, mls=trait['years'], body_mass_g=trait.get('body_mass_g'),
                available_fraction=fraction, region_eligible=fraction is not None and fraction >= .5,
                activity_label='Median peak count'))
    rows.sort(key=lambda row: (row['gene'], row['species']))
    (SUPPORT / 'species_medians.json').write_text(json.dumps(rows, indent=2, allow_nan=False))
    (SUPPORT / 'genes.txt').write_text('\n'.join(sorted(genes)) + '\n')
    with (ORIGINAL / 'all_tested.tsv').open() as f:
        old_rows = [r for r in csv.DictReader(f, delimiter='\t') if r['gene'] in genes]
    assert len(old_rows) == 71
    (SUPPORT / 'original_71_gene_statistics.json').write_text(json.dumps(old_rows, indent=2))
    for name in ('peaks500kb_without_body_mass_all_tested.tsv.gz',
                 'peaks500kb_with_body_mass_all_tested.tsv.gz'):
        shutil.copy2(SOURCE / name, SUPPORT / name)
    shutil.copy2(ROOT / 'cross_species_visual_tool/plot_regulatory_mls.R', SUPPORT / 'plot_regulatory_mls.R')
    provenance = dict(gene_list=str(ORIGINAL / 'FDR05_positive_genes.txt'),
        gene_list_sha256=sha256(ORIGINAL / 'FDR05_positive_genes.txt'), n_genes=len(genes),
        plot='6 x 6 inches, 1800 x 1800 pixels, black points, red dashed uncorrected fit',
        analysis_manifest=manifest, cache_fingerprints=cache_fingerprints,
        data_sha256=sha256(SUPPORT / 'species_medians.json'),
        fdr_family_sha256={name:sha256(SUPPORT / name) for name in
            ('peaks500kb_without_body_mass_all_tested.tsv.gz','peaks500kb_with_body_mass_all_tested.tsv.gz')})
    (SUPPORT / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    print(f'Prepared {len(genes)} genes, {len(cache_fingerprints)} species, {len(rows)} gene-species rows.')
    print('Original table columns:', list(old_rows[0]))


if __name__ == '__main__':
    main()
