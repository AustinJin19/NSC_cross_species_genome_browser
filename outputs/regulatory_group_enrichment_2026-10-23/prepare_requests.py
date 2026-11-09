"""Prepare comparable GO:BP and KEGG queries for three regulatory hit groups."""
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent.parent
REG=ROOT/'outputs/regulatory_windows_2026-10-22'
RAW=ROOT/'outputs/regulatory_windows_no_mass_2026-10-23'
REF=ROOT/'outputs/aging_longevity_reference_2026-10-22'
TABLE=ROOT/'outputs/330_genes_integrated_MLS_2026-10-23/results.json'
WINDOWS=[20,50,100,200,500,1000]

def read(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()) as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

hgnc=read(REF/'sources/hgnc_complete_set.tsv')
approved={r['symbol'].upper():r['symbol'] for r in hgnc}
hgnc_by_symbol={r['symbol']:r for r in hgnc}
aliases=defaultdict(set)
for r in hgnc:
    for a in (r['alias_symbol']+'|'+r['prev_symbol']).split('|'):
        if a:aliases[a.upper()].add(r['symbol'])
def canonical(g):
    return approved.get(g.upper()) or (next(iter(aliases[g.upper()])) if len(aliases[g.upper()])==1 else g.upper())
def identifier(g):
    c=canonical(g)
    ensembl=hgnc_by_symbol.get(c,{}).get('ensembl_gene_id','')
    return ensembl if re.fullmatch(r'ENSG\d+',ensembl) else c

tested=set()
tested_windows=defaultdict(list)
sources=[TABLE,REF/'sources/hgnc_complete_set.tsv']
for kb in WINDOWS:
    pa=REG/f'tss_{kb}kb/all_tested.tsv.gz';pb=RAW/f'tss_{kb}kb/all_tested.tsv.gz'
    a={r['gene'] for r in read(pa)};b={r['gene'] for r in read(pb)}
    assert a==b, kb
    tested.update(a)
    for g in a:tested_windows[g].append(kb)
    sources.extend([pa,pb])
background=sorted({canonical(g) for g in tested})
background_ids=sorted({identifier(g) for g in tested})
groups=defaultdict(list)
labels={'Both':'Shared','Corrected only':'Corrected only','Uncorrected only':'Uncorrected only'}
slugs={'Shared':'shared','Corrected only':'corrected_only','Uncorrected only':'uncorrected_only'}
data=json.loads(TABLE.read_text())
for row in data['overview']:
    assert row['Gene'] in tested and row['Human_symbol']==canonical(row['Gene'])
    groups[labels[row['Regulatory_significance']]].append(row)
assert {g:len(rows) for g,rows in groups.items()}=={'Shared':143,'Corrected only':73,'Uncorrected only':114}
assert len({r['Human_symbol'] for rows in groups.values() for r in rows})==330
out_genes=[]
requests=[]
for group in ['Shared','Corrected only','Uncorrected only']:
    rows=groups[group]
    genes=sorted(r['Human_symbol'] for r in rows)
    gene_ids=sorted({identifier(g) for g in genes})
    assert len(gene_ids)==len(genes)
    assert set(genes)<=set(background)
    for r in rows:
        out_genes.append(dict(Group=group,Gene=r['Gene'],Human_symbol=r['Human_symbol'],Submitted_ID=identifier(r['Gene']),
            Corrected_windows_kb=r['Corrected_windows_kb'],Uncorrected_windows_kb=r['Uncorrected_windows_kb'],
            Corrected_direction=r['Regulatory_corrected'],Uncorrected_direction=r['Regulatory_uncorrected']))
    for library,short in [('GO:BP','GO_BP'),('KEGG','KEGG')]:
        name=slugs[group]+'_'+short
        request=dict(organism='hsapiens',query=gene_ids,background=background_ids,domain_scope='custom',sources=[library],
                     significance_threshold_method='fdr',user_threshold=.05,all_results=True,
                     no_evidences=False,ordered=False,no_iea=False)
        # Retain the initial symbol-only trial, which exposed obsolete-name mappings.
        archive=OUT/'preliminary_symbol_queries'
        archive.mkdir(exist_ok=True)
        for suffix in ['request','response']:
            current=OUT/f'{name}_{suffix}.json'
            if current.exists() and not (archive/current.name).exists():
                current.rename(archive/current.name)
        (OUT/f'{name}_request.json').write_text(json.dumps(request))
        requests.append(dict(name=name,group=group,library=library,submitted_genes=len(genes)))
write(OUT/'input_genes.tsv',out_genes)
write(OUT/'background_symbol_mapping.tsv',[dict(Original_symbol=g,Human_symbol=canonical(g),
    Submitted_ID=identifier(g),HGNC_approved=canonical(g) in hgnc_by_symbol,Tested_windows_kb=';'.join(map(str,tested_windows[g]))) for g in sorted(tested)])
(OUT/'background_genes.txt').write_text('\n'.join(background_ids)+'\n')
(OUT/'manifest.json').write_text(json.dumps(dict(requests=requests,original_tested_genes=len(tested),
    canonical_background_genes=len(background),submitted_background_ids=len(background_ids),groups={g:len(rows) for g,rows in groups.items()},
    identifier_rule='Use HGNC-provided single Ensembl ID for an approved human symbol; otherwise submit symbol and leave service failures unresolved',
    background_rule='Union of genes tested in at least one of the six windows; identical tested sets in both models',
    correction='g:Profiler BH FDR separately within each of six group/library queries; extra joint BH across all returned tests will be provided',
    sources_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    api='https://biit.cs.ut.ee/gprofiler/api/gost/profile/'),indent=2))
print(json.dumps(dict(original_tested_genes=len(tested),canonical_background_genes=len(background),
    query_genes=330,groups={g:len(rows) for g,rows in groups.items()},requests=requests),indent=2))
