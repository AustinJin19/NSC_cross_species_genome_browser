import csv,json
from pathlib import Path
root=Path(__file__).resolve().parent
summary={}
for group in ['positive','negative']:
 d=json.loads((root/f'{group}_response.json').read_text());request=json.loads((root/f'{group}_request.json').read_text())
 gm=d['meta']['genes_metadata'];qm=gm['query']['query_1'];reverse={v:k for k,vs in qm['mapping'].items() for v in vs}
 rows=[]
 for r in d['result']:
  ids=[g for g,ev in zip(qm['ensgs'],r['intersections']) if ev]
  assert len(ids)==r['intersection_size']
  rows.append(dict(term_ID=r['native'],ontology=r['source'],term=r['name'],FDR_adjusted_p=r['p_value'],significant=r['significant'],overlap=r['intersection_size'],query_size=r['query_size'],term_in_background=r['term_size'],background_size=r['effective_domain_size'],fold_enrichment=(r['intersection_size']/r['query_size'])/(r['term_size']/r['effective_domain_size']),genes=';'.join(reverse.get(g,g) for g in ids)))
 rows.sort(key=lambda r:r['FDR_adjusted_p'])
 for suffix,selected in [('all_terms',rows),('significant_terms',[r for r in rows if r['significant']])]:
  with (root/f'{group}_{suffix}.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(selected)
 with (root/f'{group}_mapping.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['input_symbol','mapping_status','human_ensembl_ids'])
  for g in request['query']:
   status='mapped' if g in qm['mapping'] else 'ambiguous' if g in gm['ambiguous'] else 'unmapped'
   w.writerow([g,status,';'.join(qm['mapping'].get(g,[]))])
 summary[group]=dict(input_genes=len(request['query']),mapped_genes=len(qm['ensgs']),unmapped=len(gm['failed']),ambiguous=len(gm['ambiguous']),effective_background=rows[0]['background_size'],significant=[r for r in rows if r['significant']],version=d['meta']['version'])
(root/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
