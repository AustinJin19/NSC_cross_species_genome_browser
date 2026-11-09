"""Overlap the reference with per-window FDR-significant regulatory activity genes."""
import csv, json
from pathlib import Path
from collections import defaultdict

root=Path(__file__).resolve().parent
search=root.parent/'regulatory_windows_2026-10-22'
def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
reference=read(root/'aging_longevity_genes.tsv')
ref={r['Human gene']:r for r in reference}
hgnc=read(root/'sources/hgnc_complete_set.tsv')
approved={r['symbol'].upper():r['symbol'] for r in hgnc}
aliases=defaultdict(set)
for r in hgnc:
    for alias in (r['alias_symbol']+'|'+r['prev_symbol']).split('|'):
        if alias:aliases[alias.upper()].add(r['symbol'])

mapping=[]
def canonical(g):
    key=g.upper()
    if key in approved:return approved[key]
    if len(aliases[key])==1:return next(iter(aliases[key]))
    return key

sig=read(search/'significant_genes_all_windows.tsv')
assert len(sig)==300 and len({r['gene'] for r in sig})==216
assert all(float(r['BH_FDR_window'])<0.05 for r in sig)
details=[]
for r in sig:
    gene=canonical(r['gene'])
    if gene!=r['gene']:
        mapping.append({'Original gene':r['gene'],'Human gene':gene})
    if gene in ref:
        details.append({'window_kb':int(r['window_kb']),'original_gene':r['gene'],'human_gene':gene,
            'direction':r['direction'],'n_species':int(r['n_species']),'partial_r':float(r['partial_r']),
            'P':float(r['p_value']),'FDR_window':float(r['BH_FDR_window']),
            **{k:v for k,v in ref[gene].items() if k!='Human gene'}})

summary=[]
for window in [20,50,100,200,500,1000,'Any']:
    candidates={canonical(r['gene']) for r in sig if window=='Any' or int(r['window_kb'])==window}
    hits=candidates & set(ref)
    summary.append({'window_kb':window,'significant_genes':len(candidates),'overlap_genes':len(hits),
        'overlap_percent':100*len(hits)/len(candidates),
        **{s:sum(ref[g][s]=='Yes' for g in hits) for s in ['Experiment','GenAge','CellAge','GWAS','CellAge manipulation','CellAge expression','GWAS longevity']},
        'overlap_symbols':'; '.join(sorted(hits))})

unique=[]
for gene in sorted({r['human_gene'] for r in details}):
    rows=[r for r in details if r['human_gene']==gene]
    unique.append({'human_gene':gene,'significant_windows_kb':'; '.join(str(w) for w in sorted({r['window_kb'] for r in rows})),
        'original_gene_symbols':'; '.join(sorted({r['original_gene'] for r in rows})),
        'MLS_direction':'; '.join(sorted({r['direction'] for r in rows})),
        **{k:v for k,v in ref[gene].items() if k!='Human gene'}})
old71={g.strip() for g in (root.parent/'regulatory-500kb-mass-2026-10-19/FDR05_positive_genes.txt').read_text().splitlines() if g.strip()}
assert old71=={r['gene'] for r in sig if r['window_kb']=='500'}
assert len({r['human_gene'] for r in details if r['window_kb']==500})==summary[4]['overlap_genes']
for name,rows in [('regulatory_overlap_summary.tsv',summary),('regulatory_overlap_gene_windows.tsv',details),('regulatory_overlap_unique_genes.tsv',unique),('regulatory_500kb_overlap.tsv',[r for r in details if r['window_kb']==500])]:
    with (root/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(root/'regulatory_overlap_results.json').write_text(json.dumps({'summary':summary,'unique':unique,'symbol_changes':mapping},indent=2))
print(json.dumps({'summary':summary,'500kb_genes':[{k:r[k] for k in ['human_gene','Data resources','CellAge manipulation','CellAge expression','GWAS longevity']} for r in details if r['window_kb']==500],'symbol_changes':mapping},indent=2))
