import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

HERE=Path(__file__).resolve().parent
GO=HERE/'go'
OUT=HERE.parent.parent/'outputs/regulatory-50kb-mass-2026-10-19'
(OUT/'go').mkdir(exist_ok=True)
report=dict(groups={}, mappings=[], memberships=[])
for item in json.loads((GO/'manifest.json').read_text()):
    name=item['name']
    if not item['run']:
        report['groups'][name]=dict(summary={**item,'status':'Not run: empty gene list'}, terms=[])
        continue
    d=json.loads((GO/f'{name}_response.json').read_text())
    request=json.loads((GO/f'{name}_request.json').read_text())
    gm=d['meta']['genes_metadata']; qm=gm['query']['query_1']
    reverse=defaultdict(list)
    for symbol,ids in qm['mapping'].items():
        for gene in ids: reverse[gene].append(symbol)
    rows=[]
    for r in d['result']:
        ids=[g for g,ev in zip(qm['ensgs'],r['intersections']) if ev]
        assert len(ids)==len(set(ids))==r['intersection_size']
        N,K,n,k=r['effective_domain_size'],r['term_size'],r['query_size'],r['intersection_size']
        row=dict(list=name,go_id=r['native'],term=r['name'],go_fdr=r['p_value'],
                 raw_p=float(hypergeom.sf(k-1,N,K,n)),overlap=k,query_size=n,
                 term_in_background=K,background_size=N,fold_enrichment=(k/n)/(K/N),
                 significant=r['significant'],genes=';'.join('/'.join(reverse[g]) for g in ids))
        assert row['significant']==(row['go_fdr']<.05)
        rows.append(row)
        if row['significant']:
            report['memberships'].extend(dict(list=name,go_id=r['native'],term=r['name'],
                gene='/'.join(reverse[g]),ensembl_id=g) for g in ids)
    error=None
    if rows:
        p=np.array([r['raw_p'] for r in rows]); order=np.argsort(p); q=np.empty(len(p))
        q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        error=float(max(abs(q-np.array([r['go_fdr'] for r in rows]))))
        assert error<1e-10
    rows.sort(key=lambda r:(r['go_fdr'],r['raw_p'],r['go_id']))
    counts=defaultdict(int)
    for g in request['query']:
        if g in qm['mapping']: status='Mapped'; ids=qm['mapping'][g]
        elif g in gm['ambiguous']: status='Ambiguous; excluded'; ids=[x['gene_id'] for x in gm['ambiguous'][g]]
        elif g in gm['duplicates']: status='Duplicate mapping'; ids=[]
        else: assert g in gm['failed']; status='Unmapped; excluded'; ids=[]
        counts[status]+=1
        report['mappings'].append(dict(list=name,gene=g,status=status,human_ensembl_ids=';'.join(ids)))
    summary={**item,'status':'Completed','mapped_genes':len(qm['ensgs']),
             'mapping_counts':dict(counts),'returned_terms':len(rows),
             'effective_background':rows[0]['background_size'] if rows else None,
             'significant_terms':sum(r['significant'] for r in rows),
             'minimum_go_fdr':min((r['go_fdr'] for r in rows),default=None),
             'version':d['meta']['version'],'timestamp':d['meta']['timestamp'],
             'bh_max_absolute_error':error}
    report['groups'][name]=dict(summary=summary,terms=rows)
    fields=list(rows[0]) if rows else ['list','go_id','term','go_fdr','raw_p','overlap','query_size','term_in_background','background_size','fold_enrichment','significant','genes']
    for suffix,chosen in [('all_terms',rows),('significant_terms',[r for r in rows if r['significant']])]:
        with (OUT/'go'/f'{name}_{suffix}.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(chosen)
    print(name,json.dumps(summary))
    print('Top terms:',json.dumps([{k:r[k] for k in ['go_id','term','go_fdr','overlap','fold_enrichment']} for r in rows[:5]]))
(GO/'report_data.json').write_text(json.dumps(report))
(OUT/'go'/'summary.json').write_text(json.dumps({k:v['summary'] for k,v in report['groups'].items()},indent=2))
