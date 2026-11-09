"""Map the 17 curated worm/fly genes, preserving ambiguous orthology.

Raw DIOPT v9 endpoint responses report dataset version 9.1. The advertised
v10 endpoint returned 'Version not suppported', so no v10 results are claimed.
"""
import collections
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / 'sources/orthology'
TAXA = {'Caenorhabditis elegans': 6239, 'Drosophila melanogaster': 7227}
TARGETS = {9606: 'Human', 10090: 'Mouse'}
DATE = '2026-10-22'
data = json.loads((ROOT / 'collection.json').read_text())
genes = list(csv.DictReader((ROOT / 'curated_gene_list.tsv').open(), delimiter='\t'))
genes = [g for g in genes if g['Scientific name'] in TAXA]
ncbi = json.loads((SOURCES / 'ncbi_gene_summaries.json').read_text())['result']
biotypes = {
    g.findtext('.//Gene-track_geneid'): g.find('Entrezgene_type').attrib['value']
    for g in ET.parse(SOURCES / 'ncbi_gene_records.xml').getroot().findall('Entrezgene')
}
all_pairs = []
summary = []
overrides = {
    'daf-16': 'FOXO family mapping is not exclusive to FOXO3. FOXO3B is retained because current NCBI and HGNC classify it as protein coding; historical pseudogene labels should not be substituted.',
    'daf-2': 'One worm insulin-like receptor maps to three mammalian receptor genes. The worm experiment does not establish that inhibiting all three extends mammalian life.',
    'hlh-30': 'Retain the MITF/TFE family alternatives; the frequently used TFEB analogy is not an exclusive ortholog assignment.',
    'hsf-1': 'Multiple heat-shock-factor family matches. Mouse Hsf3 is included by direct mouse mapping rather than inferred from human capitalization.',
    'ifg-1': 'Multiple eIF4G-family targets; do not assign only EIF4G1 from the common functional shorthand.',
    'isp-1': 'UQCRFS1P1 is excluded from the recommended mapping because NCBI classifies it as a pseudogene.',
    'lipl-4': 'Many-to-many lipase-family relationship; no exclusive LIPA conversion. Mouse-specific lipase-family genes remain in the expanded mapping.',
    'pat-10': 'Human CALML6 is supported by DIOPT and the NCBI/Alliance source-gene summary. No mouse hit reaches the score cutoff; weak calmodulin-family candidates are in the audit TSV. Do not automatically substitute TNNC1/TNNC2 from troponin-C function.',
    'rsks-1': 'Both RPS6KB1 and RPS6KB2 retained; the mouse Rps6kb1 lifespan experiment is distinct from this orthology inference.',
    'skn-1': 'Retain NFE2 and NFE2L1/2/3; the common NRF2 analogy does not establish exclusive NFE2L2 orthology.',
    'xbp-1': 'The source lifespan intervention used neuronal spliced XBP-1s; orthology does not preserve the experimental tissue or isoform automatically.',
    'Atg8a': 'GABARAP-family targets pass the chosen cutoff. Lower-support LC3-family predictions are retained in the audit, not the recommended list.',
    'Pten': 'PTEN/Pten is the strongest match (13/13). TPTE/Tpte is an additional family prediction and is not interchangeable with PTEN for lifespan interpretation.',
    'Thor': 'Retain all three EIF4EBP-family genes; the beneficial fly intervention was an activating allele.',
    'chico': 'Human IRS1/2/4 and mouse Irs1/2/3/4 differ; mouse Irs3 is obtained directly from the mouse mapping.',
    'foxo': 'FOXO family mapping, including the current protein-coding FOXO3B annotation in human; no unique FOXO3 assignment.',
    'park': 'PRKN/Prkn is the strongest reciprocal match; PARK2 is a familiar older alias.',
}

for g in genes:
    symbol, organism = g['Gene'], g['Scientific name']
    taxid = TAXA[organism]
    selected = {}
    source_id = None
    versions = set()
    api_urls = {}
    for target, common in TARGETS.items():
        raw = json.loads((SOURCES / f'{taxid}_{symbol}_{target}.json').read_text())
        search = raw['search_details']
        assert search['input_species'] == taxid and search['target_species'] == target
        assert search['gene_details'][0]['symbol'].lower() == symbol.lower()
        source_id = str(search['search_gene_entrez'])
        source_stable = search['gene_details'][0]['species_specific_geneid']
        versions.add(str(search['diopt_version']))
        url = f'https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez/{taxid}/{source_id}/{target}/none'
        api_urls[target] = url
        hits = list(raw.get('results', {}).get(source_id, {}).values())
        retained = []
        for h in hits:
            assert h['species_id'] == target
            target_id = str(h['geneid'])
            annotation = ncbi.get(target_id, {})
            current_symbol = annotation.get('nomenclaturesymbol') or annotation.get('name') or h['symbol']
            biotype = biotypes.get(target_id, 'Not checked: low-support result')
            support = h['score'] >= 4 and h['confidence'] in ('high', 'moderate')
            if not support:
                decision = 'Excluded: below support threshold'
            elif biotype != 'protein-coding':
                decision = 'Excluded: not confirmed protein coding'
            elif annotation.get('status') or annotation.get('currentid'):
                decision = 'Excluded: discontinued/replaced NCBI record'
            else:
                decision = 'Retained'
            row = {
                'Source gene': symbol, 'Source species': organism,
                'Source Entrez ID': source_id, 'Source stable ID': source_stable,
                'Target species': common, 'Target gene': current_symbol,
                'Target symbol in DIOPT': h['symbol'], 'Target Entrez ID': target_id,
                'Target gene type (NCBI)': biotype,
                'DIOPT score': h['score'], 'Maximum score in response': h['max_score'],
                'DIOPT confidence': h['confidence'],
                'Best forward match': h['best_score'], 'Best reverse match': h['best_score_rev'],
                'Mapping decision': decision,
                'Supporting methods': '; '.join(h['methods']),
                'Source lifespan direction': g['Beneficial direction'],
                'Experiment records': g['Record IDs'],
                'DIOPT dataset version': str(search['diopt_version']),
                'Retrieved': DATE, 'DIOPT API URL': url,
                'Target NCBI Gene URL': f'https://www.ncbi.nlm.nih.gov/gene/{target_id}',
            }
            all_pairs.append(row)
            if decision == 'Retained':
                retained.append(row)
        selected[target] = sorted(retained, key=lambda r: (-r['DIOPT score'], r['Target gene']))
    h, m = selected[9606], selected[10090]
    summary.append({
        'Source gene': symbol,
        'Source animal': 'Worm' if taxid == 6239 else 'Fruit fly',
        'Human homolog candidates': '; '.join(r['Target gene'] for r in h) if h else 'Unresolved',
        'Mouse homolog candidates': '; '.join(r['Target gene'] for r in m) if m else 'Unresolved at score >=4',
        'Human scores': '; '.join(f"{r['Target gene']} {r['DIOPT score']}/{r['Maximum score in response']}" for r in h),
        'Mouse scores': '; '.join(f"{r['Target gene']} {r['DIOPT score']}/{r['Maximum score in response']}" for r in m),
        'Source beneficial direction': g['Beneficial direction'],
        'Experiment records': g['Record IDs'],
        'Mapping note': overrides[symbol],
        'Evidence boundary': 'Lifespan extension was tested in the source animal; mapping alone is not mammalian lifespan evidence.',
        'Mapping method': 'DIOPT score >=4; confidence high/moderate; current protein-coding NCBI record',
        'DIOPT dataset version': '; '.join(sorted(versions)),
        'Source NCBI Gene URL': f'https://www.ncbi.nlm.nih.gov/gene/{source_id}',
        'Human mapping URL': api_urls[9606],
        'Mouse mapping URL': api_urls[10090],
    })

retained = [r for r in all_pairs if r['Mapping decision'] == 'Retained']
review_pairs = [r for r in all_pairs if r['DIOPT score'] >= 4 or (r['Source gene'] == 'pat-10' and r['DIOPT confidence'] in ('high', 'moderate'))]
def tsv(name, rows):
    with (ROOT / name).open('w', encoding='utf-8', newline='') as out:
        writer = csv.DictWriter(out, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
tsv('worm_fly_to_human_mouse.tsv', summary)
tsv('ortholog_mapping_details.tsv', review_pairs)
tsv('ortholog_all_predictions_audit.tsv', all_pairs)
stats = {
    'source_genes': len(summary), 'retained_pairs': len(retained),
    'raw_prediction_pairs': len(all_pairs), 'review_pairs': len(review_pairs),
    'mapped_source_genes': {v: sum(any(r['Source gene'] == g['Gene'] and r['Target species'] == v for r in retained) for g in genes) for v in TARGETS.values()},
    'distinct_targets': {v: len({r['Target Entrez ID'] for r in retained if r['Target species'] == v}) for v in TARGETS.values()},
    'dataset_version': '9.1', 'retrieved': DATE,
    'score_filter': 'score >=4 and confidence high/moderate; current protein-coding NCBI gene',
    'scope': '17 curated C. elegans / D. melanogaster genes only; not all GenAge screening records',
}
for target in TARGETS.values():
    symbols = sorted({r['Target gene'] for r in retained if r['Target species'] == target})
    (ROOT / f'{target.lower()}_ortholog_candidates.txt').write_text('\n'.join(symbols) + '\n')
manifest = {
    **stats,
    'api_note': 'v10 endpoint was advertised but unavailable; working v9 endpoint explicitly returned dataset 9.1.',
    'sources': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SOURCES.glob('*')) if p.is_file()},
}
(ROOT / 'orthology_manifest.json').write_text(json.dumps(manifest, indent=2))
(ROOT / 'orthology.json').write_text(json.dumps({'summary': summary, 'pairs': review_pairs, 'stats': stats}, indent=2))
(ROOT / 'ORTHOLOGY_README.txt').write_text(f"""WORM / FLY TO HUMAN AND MOUSE | {DATE}

Scope: the 17 worm/fly genes in the curated primary-experiment list. The 852-record GenAge screen has NOT been mapped by this operation.

Source identifiers were matched by species and symbol to the original GenAge export and confirmed against each DIOPT query response. Human (9606) and mouse (10090) were queried independently. The source gene and beneficial perturbation direction remain attached to each prediction.

Selection: retain DIOPT score >=4 and high/moderate confidence, plus a current protein-coding NCBI Gene record. This is our stated screening threshold, not a claim of experimental validation of orthology or of a unique one-to-one relationship. Scores count supporting methods; they are NOT probabilities. The score denominator (13 in these responses) is preserved exactly. Lower-support matches and excluded pseudogenes are preserved in ortholog_all_predictions_audit.tsv.

Retrieved DIOPT dataset version 9.1 through the documented v9 endpoint. A v10 endpoint was advertised but returned 'Version not suppported'. We therefore do NOT claim v10 results. Current NCBI Gene summaries and XML records were checked for target symbols, status and biotype; the human FOXO3B record was additionally checked in HGNC (approved gene with protein product). UQCRFS1P1 was excluded as a pseudogene. Raw source responses and checksums are retained.

Results: {stats['retained_pairs']} retained source-target pairs; {stats['distinct_targets']['Human']} distinct human candidates and {stats['distinct_targets']['Mouse']} distinct mouse candidates. All 17 source genes have a retained human candidate; 16 have mouse candidates. pat-10 has no mouse candidate meeting the cutoff and remains unresolved for mouse. Its human CALML6 match agrees with the NCBI/Alliance description; troponin-C function alone does not justify converting it to TNNC1/TNNC2.

Multiple targets are intentionally retained. lipl-4 maps to a broad lipase family; daf-16 and fly foxo map to FOXO-family genes; hlh-30 maps to the MITF/TFE family. A strongest DIOPT match is not necessarily a unique functional equivalent. Direct human and mouse results differ, for example Irs3 in mouse chico mapping.

LIFESPAN EVIDENCE DOES NOT TRANSFER AUTOMATICALLY. A human/mouse homolog is a candidate for comparison with your ATAC/RNA results. Its overexpression, inhibition or knockout is not thereby proven to extend human/mouse lifespan. Source experimental direction, tissue and allele must remain distinct from the orthology prediction.

Sources:
https://www.flyrnai.org/tools/diopt/web/api
https://fgr.hms.harvard.edu/diopt-documentation
https://www.ncbi.nlm.nih.gov/gene/172127
https://rest.genenames.org/fetch/symbol/FOXO3B

Outputs: worm_fly_to_human_mouse.tsv (one row per source gene); ortholog_mapping_details.tsv (retained and flagged review candidates); ortholog_all_predictions_audit.tsv (all predictions); human_ortholog_candidates.txt and mouse_ortholog_candidates.txt (deduplicated candidates, NOT lists of proven mammalian lifespan-extension genes).
""")
assert len(summary) == 17
assert stats['mapped_source_genes'] == {'Human': 17, 'Mouse': 16}
assert not any(r['Target gene'] == 'UQCRFS1P1' for r in retained)
print(json.dumps(stats, indent=2))
for r in summary:
    print(r['Source gene'], '|', r['Human homolog candidates'], '|', r['Mouse homolog candidates'])
