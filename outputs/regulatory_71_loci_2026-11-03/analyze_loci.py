"""Map the original 71 positive regulatory hits without rerunning MLS models."""
from pathlib import Path
import csv, gzip, json, re, math, itertools, collections

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SOURCE = ROOT / 'outputs/regulatory-500kb-mass-2026-10-19'
CACHE = ROOT / 'output/regulatory_500kb_mass_2026-10-19'

def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))

hits = read_tsv(SOURCE / 'FDR05_positive.tsv')
assert len(hits) == 71 and len({h['gene'] for h in hits}) == 71
symbols = {h['gene'] for h in hits}
stats = {h['gene']: h for h in hits}
hgnc = read_tsv(ROOT / 'outputs/aging_longevity_reference_2026-10-22/sources/hgnc_complete_set.tsv')
approved = {r['symbol']: r for r in hgnc}
previous = collections.defaultdict(set)
for r in hgnc:
    for s in r['prev_symbol'].split('|'):
        previous[s].add(r['symbol'])
anchors = json.loads((ROOT / 'cross_species_visual_tool/index/Human_anchors.json').read_text())
by_symbol = collections.defaultdict(list)
for a in anchors:
    if re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)', a['chrom']):
        by_symbol[a['symbol'].upper()].append(a)

mapping = {}
for g in sorted(symbols):
    s, rule = g, 'Exact symbol'
    if g not in approved and len(previous[g]) == 1:
        s, rule = next(iter(previous[g])), 'HGNC previous symbol'
    candidates = by_symbol.get(s.upper(), [])
    mapping[g] = (s, rule, candidates)
print('Unmapped or multiple:', [(g,s,len(a)) for g,(s,r,a) in mapping.items() if len(a)!=1])

# Full gene spans are the union of transcript bodies on the same chromosome/strand.
target = {s.upper() for s,_,a in mapping.values() if len(a)==1}
bodies = {}
with gzip.open(ROOT / 'cross_species_visual_tool/index/Human_models.tsv.gz', 'rt') as f:
    for line in f:
        c, start, end, tx, s, strand, *_ = line.rstrip('\n').split('\t')
        if s.upper() not in target:
            continue
        key=(s.upper(),c,strand)
        lo,hi=bodies.get(key,(int(start),int(end)))
        bodies[key]=(min(lo,int(start)),max(hi,int(end)))

def chr_key(c):
    return int(c[3:]) if c[3:].isdigit() else {'X':23,'Y':24}.get(c[3:],99)

genes=[]
for g,(s,rule,arr) in mapping.items():
    h=stats[g]
    rec=dict(gene=g,human_symbol=s,mapping=rule,n_species=int(h['n_species']),
             partial_r=float(h['partial_r']),p_value=float(h['p_value']),BH_FDR=float(h['BH_FDR']))
    if len(arr)==1:
        a=arr[0]; lo,hi=bodies[(s.upper(),a['chrom'],a['strand'])]
        rec.update(chrom=a['chrom'],start_1based=lo+1,end_1based=hi,strand=a['strand'],
                   tss_1based=a['tss'],tss_method=a['tssMethod'],transcript=a['transcript'],
                   cytoband=approved.get(s,{}).get('location',''),
                   ucsc_url=f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={a['chrom']}%3A{max(1,a['tss']-500000)}-{a['tss']+500000}")
    else:
        rec.update(chrom='',mapping='No unique human symbol/anchor mapping',cytoband='')
    genes.append(rec)
genes.sort(key=lambda r:(chr_key(r['chrom']) if r['chrom'] else 100,r.get('tss_1based',0),r['gene']))
mapped=[r for r in genes if r['chrom']]

clusters=[]
for r in mapped:
    if not clusters or clusters[-1][-1]['chrom']!=r['chrom'] or r['tss_1based']-clusters[-1][-1]['tss_1based']>=1000000:
        clusters.append([])
    clusters[-1].append(r)
cluster_rows=[]
for n,rr in enumerate(clusters,1):
    cid=f'L{n:02d}'; lo=rr[0]['tss_1based']; hi=rr[-1]['tss_1based']
    for r in rr:
        r['cluster']=cid; r['cluster_genes']=len(rr)
    cluster_rows.append(dict(cluster=cid,chrom=rr[0]['chrom'],n_genes=len(rr),
        first_tss_1based=lo,last_tss_1based=hi,tss_span_bp=hi-lo,
        common_window_bp=max(0,1000000-(hi-lo)),
        genes='; '.join(r['gene'] for r in rr),
        cytobands='; '.join(dict.fromkeys(r['cytoband'] for r in rr)),
        min_FDR=min(r['BH_FDR'] for r in rr)))

# Recover exactly the highest-median locus (lexicographic first on tied scores).
prov={r['species']:r for r in json.loads((CACHE/'species_provenance.json').read_text())}
selected={}; species_rows=[]
for f in sorted((CACHE/'species_scores').glob('*.json.gz')):
    sid=f.name.removesuffix('.json.gz')
    chosen={}
    for r in sorted(json.load(gzip.open(f,'rt')),key=lambda r:r['locus']):
        g=r['symbol_key']; v=r['median']
        if g not in symbols or v is None or not math.isfinite(v): continue
        if g not in chosen or v>chosen[g]['median']: chosen[g]=r
    for g,r in chosen.items():
        m=re.fullmatch(r'(.+)\|(.+):(\d+)-(\d+):([+-])',r['locus'])
        assert m, r['locus']
        _,chrom,start,end,strand=m.groups(); start=int(start); end=int(end)
        tss=start if strand=='+' else end-1
        included=sid in stats[g]['species'].split(',')
        rec=dict(gene=g,species=sid,scientific_name=prov[sid]['scientific_name'],
            chrom=chrom,start_1based=start+1,end_1based=end,strand=strand,tss_1based=tss+1,
            requested_window_start0=r['region_start'],requested_window_end=r['region_end'],
            window_start0=max(0,r['region_start']),window_end=max(0,r['region_start'])+r['available_bp'],
            median_peak_count=r['median'],available_fraction=r['available_fraction'],
            included_in_original_fit=included,region_eligible=r['region_eligible'],original_locus=r['locus'])
        selected[g,sid]=rec; species_rows.append(rec)
for g in symbols:
    expected=set(stats[g]['species'].split(','))
    actual={sid for (gene,sid),r in selected.items() if gene==g and r['included_in_original_fit']}
    assert actual==expected,(g,actual^expected)

# Independent residual-correlation check ties extracted loci to the original models.
import numpy as np
max_r_error=0.0
for g in symbols:
    ids=sorted(stats[g]['species'].split(','))
    x=np.log10([prov[s]['MLS'] for s in ids])
    y=np.log10([selected[g,s]['median_peak_count'] for s in ids])
    cov=np.column_stack([np.ones(len(ids)),np.log10([prov[s]['body_mass_g'] for s in ids])])
    rx=x-cov@np.linalg.lstsq(cov,x,rcond=None)[0]
    ry=y-cov@np.linalg.lstsq(cov,y,rcond=None)[0]
    err=abs(float(np.corrcoef(rx,ry)[0,1])-float(stats[g]['partial_r']))
    max_r_error=max(max_r_error,err)
    assert err<1e-10,(g,err)
for r in species_rows:
    assert r['window_end']>r['window_start0']>=0
    assert abs((r['window_end']-r['window_start0'])/1000000-r['available_fraction'])<1e-8

pair_rows=[]
for a,b in itertools.combinations(sorted(symbols),2):
    aa=next(r for r in genes if r['gene']==a); bb=next(r for r in genes if r['gene']==b)
    human_dist=abs(aa['tss_1based']-bb['tss_1based']) if aa['chrom'] and aa['chrom']==bb['chrom'] else None
    overlaps=[]; eligible=0; same_chrom=0; details=[]; count_equal=0; count_pairs=[]
    for sid in prov:
        x=selected.get((a,sid)); y=selected.get((b,sid))
        if not x or not y or not x['included_in_original_fit'] or not y['included_in_original_fit']: continue
        eligible+=1
        count_pairs.append((x['median_peak_count'],y['median_peak_count']))
        if x['chrom']!=y['chrom']: continue
        same_chrom+=1
        ov=max(0,min(x['window_end'],y['window_end'])-max(x['window_start0'],y['window_start0']))
        if ov:
            # Denominator is shorter actually available interval, not the requested 1 Mb.
            frac=ov/min(x['window_end']-x['window_start0'],y['window_end']-y['window_start0'])
            overlaps.append(frac); details.append(sid)
            count_equal+=int(x['median_peak_count']==y['median_peak_count'])
    if (human_dist is not None and human_dist<1000000) or overlaps:
        import statistics
        pair_rows.append(dict(gene_a=a,gene_b=b,human_same_chrom=human_dist is not None,
            human_tss_distance_bp=human_dist,human_window_overlap_bp=max(0,1000000-human_dist) if human_dist is not None else 0,
            species_both_in_fit=eligible,species_same_contig=same_chrom,species_overlapping=len(overlaps),
            overlap_species_fraction=len(overlaps)/eligible if eligible else None,
            median_overlap_fraction=statistics.median(overlaps) if overlaps else None,
            overlapping_species='; '.join(details)))
pair_rows.sort(key=lambda r:(-r['species_overlapping'],r['human_tss_distance_bp'] or 1e12,r['gene_a']))

def tsv(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
tsv('human_gene_locations.tsv',genes)
tsv('human_window_clusters.tsv',cluster_rows)
tsv('original_species_loci.tsv',sorted(species_rows,key=lambda r:(r['gene'],r['species'])))
tsv('window_overlap_pairs.tsv',pair_rows)
summary=dict(total_genes=len(genes),human_mapped=len(mapped),human_unmapped=[g['gene'] for g in genes if not g['chrom']],
    human_clusters=len(clusters),multi_gene_clusters=sum(len(r)>1 for r in clusters),
    genes_in_multi_gene_clusters=sum(len(r) for r in clusters if len(r)>1),
    singleton_genes=sum(len(r)==1 for r in clusters),selected_species_loci=len(species_rows),
    species=len(prov),human_overlapping_pairs=sum(r['human_window_overlap_bp']>0 for r in pair_rows),
    maximum_original_r_reproduction_error=max_r_error)
payload=dict(summary=summary,genes=genes,clusters=cluster_rows,species_loci=species_rows,pairs=pair_rows)
(OUT/'locus_analysis.json').write_text(json.dumps(payload,indent=2))
print(json.dumps(summary,indent=2))
print('\nMulti-gene clusters:')
for r in sorted(cluster_rows,key=lambda r:-r['n_genes']):
    if r['n_genes']>1: print(r)
print('\nSIRT6 pairs:',*[p for p in pair_rows if 'SIRT6' in (p['gene_a'],p['gene_b'])],sep='\n')
