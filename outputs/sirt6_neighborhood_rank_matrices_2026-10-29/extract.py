"""Rank the 63 screenshot labels using current browser sources (read-only)."""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
WEB = ROOT/'cross_species_visual_tool/web_app'
sys.path.insert(0, str(WEB))
from server import SPECS, ROOT as INDEX_ROOT, activity_unit
from sample_qc import EXCLUDED_SAMPLES, sample_id
from calculate_gene_activity import quantify
import omics

GENES = '''SEMA6B UBXN6 SH3GL1-AS1 EBI3 SIRT6 Metazoa_SRP EEF2 ATCAY ZFR2 RAX2
PIP5K1C CACTIN-AS1 DOHH CELF5 NCLN GNA15-DT TLE5 TLE6 LRG1 MIR4746 MPND TMIGD2
ANKRD24 MAP2K2 PIAS4 SNORD37 FTLP5 MATK TJP3 CACTIN MFSD12-AS1 S1PR4 TLE2 PLIN5
CHAF1A EIF1P6 YJU2 CREB3L3 ZBTB7A DAPK3 MRPL54 TBXA2R FZR1 SMIM44 GNA15 PLIN4
SH3GL1 FSD1 MIR637 APBA3 GIPC3 SMIM24 RPL35P10 HDGFL2 STAP2 NMRK2 HMG20B NFIC
GNA11 SHD RN7SL202P MFSD12 TEKTIP1'''.split()
assert len(GENES) == len(set(g.upper() for g in GENES)) == 63
LOOKUP = {g.upper():g for g in GENES}
EXCLUDED_SPECIES = {'Human','Balaena_mysticetus'}
SPECS = sorted([s for s in SPECS if s[0] not in EXCLUDED_SPECIES], key=lambda s:s[1].lower())
ALIASES = {'Ellobius lutescens':'Bramus lutescens'}
species = [dict(id=s[0], common=s[1].replace(' (hg38)','').replace(' (mm10)',''), latin=s[2]) for s in SPECS]
for s in species:
    if s['id']=='Bramus_lutescens': s['common']='Transcaucasian mole vole'
sid_by_latin = {s['latin']:s['id'] for s in species}

def finite(x): return x is not None and math.isfinite(x)
def median(vals):
    vals = [v for v in vals if finite(v)]
    return statistics.median(vals) if vals else None
def write_tsv(name, fields, rows, compressed=False):
    path = OUT/name
    with (gzip.open(path,'wt') if compressed else path.open('w')) as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader()
        writer.writerows({k:'NA' if r.get(k) is None else r.get(k) for k in fields} for r in rows)

state=omics.query('SIRT6')
assert not state['stale'], f'RNA index needs refresh: {state["stale"]}'
marks=','.join('?' for _ in GENES)
rna_samples=[]; rna_groups=defaultdict(list); seen=set(); all_rna_species=set()
with sqlite3.connect(f'file:{omics.DB}?mode=ro',uri=True) as db:
    db.row_factory=sqlite3.Row
    all_rna_species={ALIASES.get(r[0],r[0]) for r in db.execute("SELECT DISTINCT species FROM measurements WHERE assay='rna' AND sample!='BMR_Turk_Lu'")}
    for raw in db.execute(f"SELECT * FROM measurements WHERE assay='rna' AND symbol IN ({marks}) AND sample!='BMR_Turk_Lu'", list(LOOKUP)):
        r=dict(raw); sid=sid_by_latin.get(ALIASES.get(r['species'],r['species']))
        if sid is None: continue
        key=(sid,r['symbol'],r['sample']); assert key not in seen,key
        seen.add(key)
        r.update(species_id=sid,gene=LOOKUP[r['symbol']])
        rna_samples.append(r)
        if finite(r['log']): rna_groups[(sid,r['gene'])].append(r['log'])
rna_values={(sid,g):median(v) for (sid,g),v in rna_groups.items()}

with sqlite3.connect(f'file:{WEB}/data/gene_activity.sqlite?mode=ro',uri=True) as db:
    db.row_factory=sqlite3.Row
    loci=[]
    for table in ['activity','activity_additions']:
        loci.extend(dict(r) for r in db.execute(f'SELECT * FROM {table} WHERE symbol_key IN ({marks})',list(LOOKUP)))
loci_by_species=defaultdict(list)
for r in loci:
    chrom,span,strand=r['locus'].rsplit('|',1)[1].rsplit(':',2)
    a,b=map(int,span.split('-')); tss=a if strand=='+' else b-1
    loci_by_species[r['species']].append(dict(gene=LOOKUP[r['symbol_key']],locus=r['locus'],chrom=chrom,
        region_start=tss-5000,region_end=tss+5000,scores=[],widths=[],samples=[]))
atac_values={}; selected_loci=[]; atac_sample_rows=[]; sources=[]
for sid,common,latin,color,folder,prefix in SPECS:
    paths=sorted(p for p in (INDEX_ROOT/folder).glob('*.bw')
                 if p.name.lower().startswith(prefix.lower()) and sample_id(p) not in EXCLUDED_SAMPLES)
    assert paths, sid
    assert activity_unit(paths)=='CPM_signal_bp', sid
    rows=loci_by_species[sid]
    for path in paths:
        stat=path.stat(); sources.append(dict(species=sid,path=str(path.relative_to(INDEX_ROOT.parent)) if path.is_relative_to(INDEX_ROOT.parent) else str(path.resolve()),bytes=stat.st_size,modified_ns=stat.st_mtime_ns))
        scores,represented,widths,statuses=quantify(path,rows)
        for i,r in enumerate(rows):
            value=float(scores[i]) if math.isfinite(scores[i]) else None
            r['scores'].append(value); r['widths'].append(int(widths[i])); r['samples'].append(path.stem)
            atac_sample_rows.append(dict(species=sid,gene=r['gene'],locus=r['locus'],sample=path.stem,
                activity=value,available_bp=int(widths[i]),available_fraction=int(widths[i])/10000,status=statuses[i]))
    by_gene=defaultdict(list)
    for r in rows:
        r['median']=median(r['scores'])
        fractions=[w/10000 for v,w in zip(r['scores'],r['widths']) if finite(v)]
        r['available_fraction']=min(fractions) if fractions else None
        by_gene[r['gene']].append(r)
    for gene,candidates in by_gene.items():
        # Same locus choice as the live panel: highest median, then locus ID.
        best=sorted(candidates,key=lambda r:(-r['median'] if finite(r['median']) else math.inf,r['locus']))[0]
        eligible=finite(best['median']) and best['available_fraction']>=.5
        atac_values[(sid,gene)]=best['median'] if eligible else None
        selected_loci.append(dict(species=sid,gene=gene,locus=best['locus'],locus_count=len(candidates),
            median_activity=best['median'],available_fraction=best['available_fraction'],eligible=eligible,
            rank_status='Included' if eligible else 'Less than 50% available' if finite(best['median']) else 'No measurement',
            samples=len([v for v in best['scores'] if finite(v)])))
    print(f'{sid}: {len(paths)} samples; {len(by_gene)} screenshot labels found',flush=True)

def matrix(values):
    ranks={}
    for sp in species:
        col={g:values.get((sp['id'],g)) for g in GENES}
        valid={g:v for g,v in col.items() if finite(v)}
        for g,v in valid.items():
            ranks[(sp['id'],g)]=1+sum(x>v for x in valid.values())+(sum(x==v for x in valid.values())-1)/2
    result=[]
    for g in GENES:
        rank_vals=[ranks.get((s['id'],g)) for s in species]
        present=[v for v in rank_vals if finite(v)]
        result.append(dict(gene=g,mean_rank=statistics.mean(present) if present else None,n_species=len(present),
                           ranks=rank_vals,values=[values.get((s['id'],g)) for s in species]))
    result.sort(key=lambda r:(r['mean_rank'] if r['mean_rank'] is not None else math.inf,r['gene'].upper()))
    return result

rna=matrix(rna_values); atac=matrix(atac_values)
metadata=dict(genes=GENES,species=species,excluded_species=sorted(EXCLUDED_SPECIES),
    rank_definition='Descending within each species among measured screenshot labels. Rank 1 is highest. Exact ties receive average ranks. NA is unavailable, never zero.',
    row_order='Ascending mean of available within-species ranks; species count shown. This ordering is descriptive and depends on data coverage.',
    expression='Median across samples of the source liver log.tpm values, without an additional logarithm.',
    activity='Median integrated CPM signal over TSS ±5 kb (10 kb total). Highest-median locus selected before the 50% availability filter. No distance weighting, body-mass or gene-length correction.',
    matching='Case-insensitive exact gene symbols as in the browser; no gene-alias substitution or orthology verification.',
    species_aliases=ALIASES,
    notes=['RNA sample BMR_Turk_Lu excluded because it is lung.',
           'African spiny mouse RNA species field is Acomys cahirinus, but the source common field says Acomys russatus. The source species field is used; identity remains unresolved.',
           'Nine browser species lack matching liver RNA in the supplied index, including macaque; their RNA cells are NA.',
           'Matrices include species without a usable MLS because this is an abundance ranking, not an MLS fit.'],
    rna_sources=state['sources'],atac_sources=sources,
    rna_species_with_source=sum(s['latin'] in all_rna_species for s in species),
    rna_genes_with_values=sum(r['n_species']>0 for r in rna),
    atac_genes_with_values=sum(r['n_species']>0 for r in atac))
metadata['notes'][2]=metadata['notes'][2].replace('Nine',str(len(species)-metadata['rna_species_with_source']))
payload={'metadata':metadata,'rna':rna,'atac':atac}
(OUT/'matrix_data.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
write_tsv('RNA_sample_values.tsv.gz',['species_id','gene','sample','raw','log','workbook','source_row','common'],rna_samples,True)
write_tsv('ATAC_sample_values.tsv.gz',['species','gene','locus','sample','activity','available_bp','available_fraction','status'],atac_sample_rows,True)
write_tsv('ATAC_selected_loci.tsv',['species','gene','locus','locus_count','median_activity','available_fraction','eligible','rank_status','samples'],selected_loci)
for label,rows in [('RNA_expression',rna),('TSS_5kb_activity',atac)]:
    fields=['Gene','Mean rank','Species measured']+[s['common'] for s in species]
    write_tsv(label+'_rank_matrix.tsv',fields,[dict(zip(fields,[r['gene'],r['mean_rank'],r['n_species'],*r['ranks']])) for r in rows])
print(json.dumps({k:metadata[k] for k in ['rna_species_with_source','rna_genes_with_values','atac_genes_with_values']},indent=2))
print('RNA top 5:',[r['gene'] for r in rna[:5]])
print('ATAC top 5:',[r['gene'] for r in atac[:5]])
