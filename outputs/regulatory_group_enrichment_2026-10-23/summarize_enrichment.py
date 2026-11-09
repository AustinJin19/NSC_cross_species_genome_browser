"""Summarize g:Profiler enrichment and independently verify each BH correction."""
import csv
import gzip
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom,false_discovery_control

OUT=Path(__file__).resolve().parent
manifest=json.loads((OUT/'manifest.json').read_text())
GROUPS=['Shared','Corrected only','Uncorrected only']
LIBRARIES=['GO:BP','KEGG']

def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows,fields=None):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(rows)

genes=read(OUT/'input_genes.tsv')
gene_by_id={r['Submitted_ID']:r for r in genes}
assert len(gene_by_id)==330
terms,summary,mapping,validation=[],[],[],[]
source_hashes={}
versions=set()
term_genes={}
for item in manifest['requests']:
    name=item['name'];group=item['group'];library=item['library']
    request=json.loads((OUT/f'{name}_request.json').read_text())
    p=OUT/f'{name}_response.json';data=json.loads(p.read_text())
    source_hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    meta=data['meta'];versions.add(meta['version'])
    gm=meta['genes_metadata'];qm=gm['query']['query_1']
    reverse=defaultdict(list)
    for submitted,ids in qm['mapping'].items():
        assert len(ids)==1,(submitted,ids)
        for gene in ids:reverse[gene].append(gene_by_id[submitted]['Human_symbol'])
    assert len(qm['ensgs'])==len(set(qm['ensgs']))
    rows=[]
    for r in data['result']:
        assert r['source']==library
        ids=[gene for gene,ev in zip(qm['ensgs'],r['intersections']) if ev]
        assert len(ids)==len(set(ids))==r['intersection_size']
        human=sorted(s for gene in ids for s in reverse[gene])
        assert len(human)==len(ids)
        N,K,n,k=r['effective_domain_size'],r['term_size'],r['query_size'],r['intersection_size']
        raw=float(hypergeom.sf(k-1,N,K,n))
        row=dict(Group=group,Library=library,Term_ID=r['native'],Term=r['name'],
                 Raw_P=raw,BH_FDR=float(r['p_value']),Overlap=k,Query_size=n,
                 Term_background=K,Background_size=N,Fold_enrichment=(k/n)/(K/N),
                 Gene_ratio=k/n,Significant_FDR05=bool(r['p_value']<.05),Genes='; '.join(human),
                 Source_URL='https://www.ebi.ac.uk/QuickGO/term/'+r['native'] if library=='GO:BP' else 'https://www.kegg.jp/pathway/hsa'+r['native'].split(':')[-1])
        assert row['Significant_FDR05']==r['significant']
        rows.append(row)
        term_genes[group,r['native']]=set(human)
    q=false_discovery_control([r['Raw_P'] for r in rows])
    max_error=float(np.max(np.abs(q-np.array([r['BH_FDR'] for r in rows]))))
    assert max_error<1e-10,(name,max_error)
    assert len({r['Background_size'] for r in rows})==1
    assert len({r['Query_size'] for r in rows})==1
    counts=Counter()
    for submitted in request['query']:
        if submitted in qm['mapping']:
            status='Mapped';ids=qm['mapping'][submitted]
        elif submitted in gm['ambiguous']:
            status='Ambiguous; excluded';ids=[r['gene_id'] for r in gm['ambiguous'][submitted]]
        elif submitted in gm['duplicates']:
            status='Duplicate mapping';ids=[]
        else:
            assert submitted in gm['failed'],submitted
            status='Unmapped; excluded';ids=[]
        counts[status]+=1
        mapping.append(dict(Group=group,Library=library,Gene=gene_by_id[submitted]['Gene'],
            Human_symbol=gene_by_id[submitted]['Human_symbol'],Submitted_ID=submitted,
            Status=status,Mapped_Ensembl_IDs='; '.join(ids)))
    rows.sort(key=lambda r:(r['BH_FDR'],r['Raw_P'],r['Term_ID']))
    s=dict(Group=group,Library=library,Submitted_genes=item['submitted_genes'],
           Mapped_genes=len(qm['ensgs']),Tested_terms=len(rows),Background_size=rows[0]['Background_size'],
           Significant_terms=sum(r['Significant_FDR05'] for r in rows),
           Minimum_FDR=min(r['BH_FDR'] for r in rows),Mapping_counts=dict(counts))
    summary.append(s);terms.extend(rows)
    validation.append(dict(Query=name,Maximum_BH_difference=max_error,Version=meta['version'],Timestamp=meta['timestamp']))
    print(json.dumps(s),flush=True)
    print(json.dumps([{k:r[k] for k in ['Term_ID','Term','BH_FDR','Overlap','Fold_enrichment','Genes']} for r in rows[:8]]),flush=True)

assert len(versions)==1
assert len({s['Background_size'] for s in summary})==1,'Background changed between groups/libraries'
for r,q in zip(terms,false_discovery_control([r['Raw_P'] for r in terms])):
    r['Joint_BH_FDR_six_queries']=float(q)
for s in summary:
    s['Joint_significant_terms']=sum(r['Joint_BH_FDR_six_queries']<.05 for r in terms if r['Group']==s['Group'] and r['Library']==s['Library'])
terms.sort(key=lambda r:(LIBRARIES.index(r['Library']),GROUPS.index(r['Group']),r['BH_FDR'],r['Raw_P'],r['Term_ID']))
significant=[r for r in terms if r['Significant_FDR05']]

# Select representative GO terms for the figure; retain every term in the tables.
# Within each group, skip near-identical gene sets (Jaccard >=0.8), max six terms.
selected=[]
for library in LIBRARIES:
    for group in GROUPS:
        chosen=[]
        for r in significant:
            if r['Library']!=library or r['Group']!=group:continue
            g=term_genes[group,r['Term_ID']]
            if any(len(g&old)/len(g|old)>=.8 for old in chosen):continue
            chosen.append(g)
            key=(library,r['Term_ID'])
            if key not in selected:selected.append(key)
            if len(chosen)>=6:break
lookup={(r['Group'],r['Library'],r['Term_ID']):r for r in terms}
comparison=[]
for library,term in selected:
    reference=next(r for r in terms if r['Library']==library and r['Term_ID']==term)
    for group in GROUPS:
        r=lookup.get((group,library,term))
        comparison.append(r if r else dict(Group=group,Library=library,Term_ID=term,Term=reference['Term'],
            Raw_P=None,BH_FDR=None,Overlap=0,Query_size=next(s['Mapped_genes'] for s in summary if s['Group']==group and s['Library']==library),
            Term_background=reference['Term_background'],Background_size=reference['Background_size'],Fold_enrichment=0,
            Gene_ratio=0,Significant_FDR05=False,Genes='',Source_URL=reference['Source_URL'],Joint_BH_FDR_six_queries=None))

for r in genes:
    statuses={x['Status'] for x in mapping if x['Gene']==r['Gene']}
    assert len(statuses)==1
    r['Mapping_status']=next(iter(statuses))
    r['Mapped_Ensembl_IDs']=next(x['Mapped_Ensembl_IDs'] for x in mapping if x['Gene']==r['Gene'])
notes=[
    dict(Item='Groups',Details='143 shared genes are significant in at least one regulatory window in both models; 73 corrected-only and 114 uncorrected-only genes are significant only in that model. Groups pool all six windows and both correlation directions. Overlap need not occur at the same window.'),
    dict(Item='Regulatory selection',Details='Gene-level BH FDR <0.05 separately within each ±20, ±50, ±100, ±200, ±500 or ±1000 kb window. log10(species-median peak count) versus log10(MLS); corrected model adds log10(body mass). Human, macaque, bowhead whale excluded; >=50% region availability and >=5 species required.'),
    dict(Item='Background',Details=f"Same background for all six enrichment queries: {manifest['original_tested_genes']:,} original genes tested in at least one window, resolving to {manifest['canonical_background_genes']:,} canonical symbols and {manifest['submitted_background_ids']:,} submitted IDs. Effective mapped background: {summary[0]['Background_size']:,}. Tested-gene sets are identical between the two correlation models at each window."),
    dict(Item='Identifier handling',Details='Use unique HGNC-approved human symbols/previous names, then HGNC-provided Ensembl IDs where available. This rescues renamed genes and avoids readthrough-symbol ambiguity. MS4A20 and ZFP865 have no usable human mapping and remain excluded. Human annotation is a functional proxy, not verified one-to-one orthology.'),
    dict(Item='Enrichment test',Details='Unordered one-sided hypergeometric overrepresentation with g:Profiler, Homo sapiens, custom background over all known genes, GO biological processes and KEGG queried separately. Electronic GO annotations included. Raw P reconstructed independently and overlaps verified from returned gene evidence.'),
    dict(Item='Primary FDR',Details='Benjamini-Hochberg FDR <0.05 within each group/library query. All returned terms, including nonsignificant terms, are retained. Service-adjusted P values were independently verified against BH from reconstructed raw P. These are six separate correction families.'),
    dict(Item='Joint FDR sensitivity',Details='A secondary column pools all returned term tests across the three groups and two libraries for a single BH correction. It is a different correction family and may increase or decrease individual adjusted P values.'),
    dict(Item='Fold enrichment',Details='(overlap genes / mapped query genes) divided by (background genes annotated to term / effective background genes). This measures overrepresentation, not pathway activity.'),
    dict(Item='Comparison plot',Details='Shows up to six representative FDR-significant terms per group/library, skipping near-identical within-group overlap gene sets with Jaccard similarity >=0.8. All significant terms are retained in the tables. Dot size is overlap count; color is -log10(FDR); filled dots have FDR <0.05 and open circles do not. A dash means no returned overlap term.'),
    dict(Item='Group interpretation',Details='Different significance across groups is not a direct test of differences between groups. Lists have different sizes and statistical power. Synaptic/neural term names are human functional annotations of genes selected from liver ATAC; they do not establish neuronal activity in liver.'),
    dict(Item='Dependence',Details='GO terms are nested and share genes. Regulatory windows overlap and nearby genes can share peaks, so gene-list and enrichment tests are not fully independent. These results do not establish causal longevity pathways.'),
    dict(Item='Annotation version',Details='g:Profiler '+next(iter(versions))+'; responses retrieved 23 October 2026. Original requests, raw responses, identifier mapping, full tested terms and checksums are saved.'),
    dict(Item='g:Profiler documentation',Details='https://biit.cs.ut.ee/gprofiler/page/apis'),
]
report=dict(summary=summary,significant=significant,comparison=comparison,input_genes=genes,mapping=mapping,
            notes=notes,version=next(iter(versions)),validation=validation,selected_terms=[list(x) for x in selected],
            background=manifest['canonical_background_genes'],effective_background=summary[0]['Background_size'])
write(OUT/'all_terms.tsv.gz',terms)
write(OUT/'significant_terms.tsv',significant,list(terms[0]))
write(OUT/'comparison_terms.tsv',comparison,list(terms[0]))
write(OUT/'enrichment_summary.tsv',[{k:v for k,v in r.items() if k!='Mapping_counts'} for r in summary])
write(OUT/'gene_mapping.tsv',mapping)
write(OUT/'input_genes_with_mapping.tsv',genes)
write(OUT/'methods.tsv',notes)
(OUT/'report_data.json').write_text(json.dumps(report,indent=2))
(OUT/'validation.json').write_text(json.dumps(dict(queries=validation,response_sha256=source_hashes,
    all_tested_terms=len(terms),significant_terms=len(significant),unique_significant_terms=len({r['Term_ID'] for r in significant}),
    mapped_queries={s['Group']:s['Mapped_genes'] for s in summary},failed=[r for r in mapping if r['Status']!='Mapped']),indent=2))
print('Significant terms:',len(significant),'Unique terms:',len({r['Term_ID'] for r in significant}),'Figure terms:',len(selected))
