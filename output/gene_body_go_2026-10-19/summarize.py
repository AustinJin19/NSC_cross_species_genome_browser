"""Summarize the two saved g:Profiler responses and independently verify BH FDR."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT = ROOT / 'outputs/go-gene-body-2026-10-19'
OUT.mkdir(parents=True, exist_ok=True)
SELECTED = ['GO:0007399', 'GO:0030030', 'GO:0007264', 'GO:0007409',
            'GO:0050808', 'GO:0030029', 'GO:0006325', 'GO:0007155']
report = {'groups': {}, 'mapping': [], 'representative_ids': SELECTED}


def tsv(path, records, columns):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        w.writeheader()
        w.writerows(records)


for group, expected in [('positive', 836), ('negative', 1936)]:
    request_file = HERE / f'{group}_request.json'
    response_file = HERE / f'{group}_response.json'
    request = json.loads(request_file.read_text())
    response = json.loads(response_file.read_text())
    assert len(request['query']) == expected
    assert len(set(request['background'])) == 20200
    assert set(request['query']) <= set(request['background'])
    gm = response['meta']['genes_metadata']
    qm = gm['query']['query_1']
    reverse = defaultdict(list)
    for symbol, ids in qm['mapping'].items():
        for gene_id in ids:
            reverse[gene_id].append(symbol)
    rows = []
    for r in response['result']:
        ids = [g for g, evidence in zip(qm['ensgs'], r['intersections']) if evidence]
        assert len(ids) == len(set(ids)) == r['intersection_size']
        assert all(g in reverse for g in ids)
        N, K, n, k = (r['effective_domain_size'], r['term_size'], r['query_size'], r['intersection_size'])
        assert 0 < k <= min(n, K) and n <= N and K <= N
        assert r['source'] == 'GO:BP'
        assert r['significant'] == (r['p_value'] < 0.05)
        rows.append(dict(go_id=r['native'], term=r['name'], fdr=r['p_value'],
                         raw_p=float(hypergeom.sf(k - 1, N, K, n)), overlap=k,
                         query_size=n, term_in_background=K, background_size=N,
                         fold_enrichment=(k / n) / (K / N),
                         significant=r['significant'],
                         genes=';'.join('/'.join(reverse[g]) for g in ids),
                         ensembl_ids=';'.join(ids)))
    p = np.array([r['raw_p'] for r in rows])
    order = np.argsort(p)
    q = np.empty(len(p))
    q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    error = float(np.max(np.abs(q - np.array([r['fdr'] for r in rows]))))
    assert error < 1e-10
    rows.sort(key=lambda r: (r['fdr'], r['raw_p'], r['go_id']))
    significant = [r for r in rows if r['significant']]
    counts = defaultdict(int)
    for symbol in request['query']:
        if symbol in qm['mapping']:
            status = 'Mapped'
            candidates = qm['mapping'][symbol]
        elif symbol in gm['ambiguous']:
            status = 'Ambiguous; excluded'
            candidates = [x['gene_id'] for x in gm['ambiguous'][symbol]]
        elif symbol in gm['duplicates']:
            status = 'Duplicate mapping; not counted again'
            candidates = []
        else:
            assert symbol in gm['failed']
            status = 'Unmapped; excluded'
            candidates = []
        counts[status] += 1
        report['mapping'].append(dict(group=group, symbol=symbol, status=status, human_ensembl_ids=';'.join(candidates)))
    summary = dict(input_genes=expected, mapped_genes=len(qm['ensgs']),
                   mapping_counts=dict(counts), submitted_background=20200,
                   effective_background=rows[0]['background_size'],
                   returned_terms=len(rows), significant_terms=len(significant),
                   minimum_fdr=min(r['fdr'] for r in rows), version=response['meta']['version'],
                   api_timestamp=response['meta']['timestamp'], bh_max_absolute_error=error,
                   request_sha256=hashlib.sha256(request_file.read_bytes()).hexdigest(),
                   response_sha256=hashlib.sha256(response_file.read_bytes()).hexdigest())
    report['groups'][group] = dict(summary=summary, terms=rows)
    for suffix, selected in [('all_terms', rows), ('significant_terms', significant)]:
        tsv(OUT / f'{group}_{suffix}.tsv', selected, list(rows[0]))
    print(group, json.dumps(summary))

assert not (set(json.loads((HERE / 'positive_response.json').read_text())['meta']['genes_metadata']['query']['query_1']['ensgs']) &
            set(json.loads((HERE / 'negative_response.json').read_text())['meta']['genes_metadata']['query']['query_1']['ensgs']))
assert report['groups']['positive']['summary']['significant_terms'] == 0
assert report['groups']['negative']['summary']['significant_terms'] == 212
assert all(next(r for r in report['groups']['negative']['terms'] if r['go_id'] == gid)['significant'] for gid in SELECTED)
tsv(OUT / 'gene_mapping.tsv', report['mapping'], ['group', 'symbol', 'status', 'human_ensembl_ids'])
(HERE / 'report_data.json').write_text(json.dumps(report))
(OUT / 'summary.json').write_text(json.dumps({g: v['summary'] for g, v in report['groups'].items()}, indent=2))
(OUT / 'METHODS.txt').write_text('''GO Biological Process overrepresentation analysis, 19 October 2026.
Positive input: 836 genes. Negative input: 1,936 genes. Both selected at nominal association p < 0.05.
Upstream model: log10(species-median gene body + 2 kb promoter ATAC activity) ~ log10(MLS) + log10(gene length) + log10(body mass).
Human, macaque and bowhead excluded. Lower-bound-only MLS excluded. At least 50% of requested genomic span available. Highest-median locus selected before availability filtering.
No input gene passes association-level BH FDR < 0.05. This enrichment is exploratory.
Background: all 20,200 genes tested in that same correlation analysis; effective mapped background 19,350.
Annotations: human symbol-based functional proxy, not verified cross-species orthology. g:Profiler version e114_eg62_p19_9dc92ac8.
API: https://biit.cs.ut.ee/gprofiler/api/gost/profile/
Organism hsapiens; sources GO:BP; domain_scope custom; unordered queries; electronic annotations included; overrepresentation only.
Hypergeometric test. Benjamini-Hochberg FDR separately for each directional GO query, as returned by g:Profiler. No additional correction across the two queries.
All returned term p-values and BH adjustments independently reproduced from reported contingency counts to absolute error < 1e-10.
Positive: 801 mapped, 26 unmapped, 9 ambiguous; 0 terms at FDR < 0.05 among 5,986 returned terms.
Negative: 1,894 mapped, 28 unmapped, 13 ambiguous, 1 duplicate mapping; 212 terms at FDR < 0.05 among 8,326 returned terms.
Ambiguous and unmapped symbols excluded; duplicate human gene mapping counted once. Mapping diagnostics retained.
Fold enrichment = (overlap / mapped query size) / (term count in background / effective background size).
Raw p-values, counts and FDR are static exported analysis results; rerun the saved API requests and summarization script to refresh.
Positive and negative refer to the adjusted activity-lifespan association, not gene expression or pathway activation.
GO parent/child terms overlap and 212 terms do not imply 212 independent pathways. Selected plotted terms illustrate themes, not a formal redundancy reduction.
No phylogenetic correction was performed in the upstream association model. Annotation and cross-species measurement differences may influence results.
GO guidance: https://geneontology.org/docs/go-enrichment-analysis/
g:Profiler documentation: https://biit.cs.ut.ee/gprofiler/page/apis
''')
