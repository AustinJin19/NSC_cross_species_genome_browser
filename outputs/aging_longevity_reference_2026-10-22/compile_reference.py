"""Compile source-labelled human candidates; retain original experimental organisms."""
import csv, gzip, hashlib, io, json, re, zipfile
from collections import defaultdict, Counter
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'sources'
OLD = ROOT.parent / 'experimental_longevity_genes_2026-10-22'
DATE = '2026-10-22'

def read_tsv(path):
    with open(path, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def archived(archive, name, delim=','):
    with zipfile.ZipFile(SRC / archive) as z:
        return list(csv.DictReader(io.StringIO(z.read(name).decode('utf-8-sig')), delimiter=delim))

def write_tsv(name, rows):
    if not rows:
        return
    with (ROOT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

def joined(vals):
    return '; '.join(sorted({str(v) for v in vals if v not in ('', None)}))

hgnc = read_tsv(SRC / 'hgnc_complete_set.tsv')
by_id = {r['entrez_id']: r for r in hgnc if r['entrez_id']}
by_symbol = {r['symbol']: r for r in hgnc}
by_ensembl = {r['ensembl_gene_id']: r for r in hgnc if r['ensembl_gene_id']}
aliases = defaultdict(set)
for r in hgnc:
    for alias in (r['alias_symbol'] + '|' + r['prev_symbol']).split('|'):
        if alias: aliases[alias].add(r['symbol'])

unresolved = []
normalization = []
def resolve(symbol='', entrez='', source=''):
    symbol, entrez = str(symbol).strip(), str(entrez).strip()
    r, method = None, ''
    if entrez and entrez in by_id: r, method = by_id[entrez], 'Entrez ID'
    elif symbol in by_symbol: r, method = by_symbol[symbol], 'Approved HGNC symbol'
    elif symbol in by_ensembl: r, method = by_ensembl[symbol], 'Ensembl ID'
    elif len(aliases.get(symbol, [])) == 1:
        r, method = by_symbol[next(iter(aliases[symbol]))], 'Unique HGNC previous/alias symbol'
    if r is not None:
        if symbol != r['symbol']:
            normalization.append({'Source':source,'Original symbol':symbol,'Original Entrez':entrez,
                                  'Human symbol':r['symbol'],'HGNC ID':r['hgnc_id'],'Rule':method})
        return r['symbol']
    if symbol.upper() not in ('', '-', 'NR', 'INTERGENIC', 'NA', 'N/A', 'NONE'):
        unresolved.append({'Source':source,'Original symbol':symbol,'Original Entrez':entrez,
                           'Reason':'No unique current HGNC match; not forced into the human list'})
    return None

# HCOP mappings are a screening layer, not evidence of a conserved lifespan effect.
species_files = {
    'mouse':'Mus musculus', 'c.elegans':'Caenorhabditis elegans',
    'fruitfly':'Drosophila melanogaster', 's.cerevisiae':'Saccharomyces cerevisiae',
    's.pombe':'Schizosaccharomyces pombe', 'zebrafish':'Danio rerio'}
hcop_id, hcop_symbol = defaultdict(list), defaultdict(list)
for tag, species in species_files.items():
    with gzip.open(SRC / f'human_{tag}_hcop.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            sources = sorted(set(re.split(r'[,;|]', r['support'])))
            if len(sources) < 3: continue
            human = resolve(r['human_symbol'], r['human_entrez_gene'], 'HCOP')
            if not human: continue
            # All 15-column files place the other species Entrez ID at column 8 and symbol at column 12.
            keys = list(r)
            key_id, key_symbol = r[keys[7]], r[keys[11]]
            item = (human, 'HCOP >=3 databases', '; '.join(sources), len(sources))
            if key_id != '-': hcop_id[(species,key_id)].append(item)
            if key_symbol != '-': hcop_symbol[(species,key_symbol.upper())].append(item)

official = defaultdict(list)
for r in read_tsv(SRC / 'genage_human_homologs.csv'):
    human = resolve(r['Symbol'], r['Entrez ID'], 'GenAge homolog export')
    if human:
        official[(r['Model Organism'],r['Model Organism Entrez ID'])].append(
            (human, 'GenAge supplied homolog', 'HAGR export; support score not supplied', ''))

def model_maps(species, entrez, symbol, use_hagr=True):
    maps = list(hcop_id.get((species,entrez), []))
    if not maps: maps += hcop_symbol.get((species,symbol.upper()), [])
    if use_hagr: maps += official.get((species,entrez), [])
    # Several HCOP rows can describe the same gene pair.
    return sorted(set(maps), key=lambda x: (x[0],x[1]))

evidence, mapping_audit = [], []
def add(human, source, subtype, record, organism, original, phenotype, direction,
        mapping, detail, url, pmid='', p=None, study='', variant=''):
    if not human: return
    evidence.append({'Human gene':human,'Resource':source,'Evidence type':subtype,
        'Source record':record,'Tested organism':organism,'Original gene':original,
        'Phenotype':phenotype,'Direction in source':direction,'Mapping':mapping,
        'Details':detail,'P value':p,'Study accession':study,'Variant':variant,
        'PMID':pmid,'Source URL':url})

def record_maps(record, species, gene, entrez, maps):
    for human, method, support, count in maps:
        mapping_audit.append({'Source record':record,'Source species':species,'Source gene':gene,
            'Source Entrez':entrez,'Human gene':human,'Mapping method':method,
            'Supporting databases':support,'Support count':count})

# Previously curated animal lifespan experiments. Worm/fly DIOPT decisions are preserved.
collection = json.loads((OLD / 'collection.json').read_text())
diopt = defaultdict(list)
for r in json.loads((OLD / 'orthology.json').read_text())['pairs']:
    if r['Target species'] == 'Human' and r['Mapping decision'] == 'Retained':
        human = resolve(r['Target gene'],r['Target Entrez ID'],'Curated DIOPT')
        if human:
            diopt[(r['Source species'],r['Source gene'])].append((human,'DIOPT 9.1 retained >=4',r['Supporting methods'],r['DIOPT score']))
experiments = []
for r in collection['experiments']:
    species, gene = r['Scientific name'], r['Gene in tested organism']
    maps = diopt.get((species,gene),[]) if species != 'Mus musculus' else model_maps(species,'',gene,False)
    record_maps(r['Record ID'],species,gene,'',maps)
    experiments.append({'Human candidates':joined(m[0] for m in maps), 'Mapping status':'Mapped' if maps else 'Not mapped in this collection', **r})
    for human in sorted({m[0] for m in maps}):
        add(human,'Experiment','Animal lifespan extension',r['Record ID'],species,gene,
            'Lifespan',r['Beneficial perturbation direction'],joined(m[1] for m in maps if m[0]==human),
            f"{r['Intervention']}. {r['Tissue and timing']}. {r['Reported lifespan result']} {r['Interpretation / limits']}",
            r['Primary source URL'])

human_genage = archived('genage_human.zip','genage_human.csv')
for r in human_genage:
    human = resolve(r['symbol'],r['entrez gene id'],'GenAge human')
    add(human,'GenAge','Human candidate','GAH'+r['GenAge ID'],'Homo sapiens',r['symbol'],
        'Aging','Not assigned','Human Entrez/HGNC', 'GenAge inclusion categories: '+r['why'],
        'https://genomics.senescence.info/genes/entry.php?hgnc='+quote(r['symbol']))

models = []
for r in archived('genage_models.zip','genage_models.csv'):
    rec = 'GAM'+r['GenAge ID']; species=r['organism']; gene=r['symbol']; entrez=r['entrez gene id']
    maps = model_maps(species,entrez,gene)
    record_maps(rec,species,gene,entrez,maps)
    url = 'https://genomics.senescence.info/genes/details.php?id='+r['GenAge ID']
    fungal = species in ('Saccharomyces cerevisiae','Schizosaccharomyces pombe','Podospora anserina')
    models.append({'Record':rec,'Original gene':gene,'Species':species,'Entrez ID':entrez,
        'Gene name':r['name'],'Lifespan effect':r['lifespan effect'],
        'Longevity influence':r['longevity influence'],
        'Human candidates':joined(m[0] for m in maps),'Mapping status':'Mapped' if maps else 'No retained human mapping',
        'Source URL':url})
    for human in sorted({m[0] for m in maps}):
        add(human,'GenAge','Fungal model homolog' if fungal else 'Animal model homolog',rec,species,gene,
            'Aging / lifespan / fitness',r['longevity influence'],joined(m[1] for m in maps if m[0]==human),
            'Source lifespan effect: '+r['lifespan effect']+'. Evidence was measured in the source organism; a human homolog is a candidate.',url)

cellage = archived('cellage.zip','cellage3.tsv','\t')
for i,r in enumerate(cellage,1):
    human=resolve(r['Gene symbol'],r['Entrez ID'],'CellAge manipulation')
    add(human,'CellAge','Senescence manipulation',f'CA{i:04d}','Homo sapiens',r['Gene symbol'],
        'Cellular senescence',r['Senescence Effect'],'Human Entrez/HGNC',
        'Senescence type: '+r['Type of senescence']+'; cancer cell: '+r['Cancer Cell']+'. Senescence effect is the database gene annotation, not a lifespan effect or a specific KO/OE instruction.',
        'https://pubmed.ncbi.nlm.nih.gov/'+r['Reference'].split(',')[0].strip()+'/',r['Reference'])
signatures = archived('cellage_signatures.zip','signatures1.csv',';')
for i,r in enumerate(signatures,1):
    human=resolve(r['gene_symbol'],r['entrez_id'],'CellAge expression')
    direction='Upregulated in senescence' if r['ovevrexp']=='1' else 'Downregulated in senescence' if r['underexp']=='1' else 'See source'
    add(human,'CellAge','Senescence expression',f'CAX{i:04d}','Homo sapiens',r['gene_symbol'],
        'Replicative senescence',direction,'Human Entrez/HGNC',
        'Expression meta-analysis; observational signature, not manipulation evidence. Source p value retained without reinterpretation.',
        'https://genomics.senescence.info/cells/signatures.php',p=float(r['p_value']))

# Classify a deliberately scoped set of longevity, healthspan and age-clock traits.
def trait_category(r):
    t=r['DISEASE/TRAIT'].lower(); m={x.strip().lower() for x in r['MAPPED_TRAIT'].split(',')}
    if 'reproductive lifespan' in t: return '', 'Reproductive span is not organism lifespan'
    if 'parental longevity' in m or 'life span determination trait' in m or 'age at death' in m:
        if 'healthspan' in m: return 'Combined aging traits','Joint longevity/healthspan analysis'
        return 'Longevity / lifespan','Direct or parental longevity, lifespan or age at death'
    if 'healthspan' in m or t in ('healthy aging','intrinsic capacity'):
        return 'Healthspan / healthy aging','Healthspan or global healthy-aging phenotype'
    if 'aging' in m and t.startswith(('biological age','epigenetic age','extrinsic epigenetic','intrinsic epigenetic','dna methylation grimage','dna methylation hannum','dna methylation phenoage','phenotypic age','phenoageaccel','physiological aging','aging (eyeage')):
        return 'Biological / epigenetic aging','Age estimate, age acceleration or aging rate; includes stated subgroup context'
    if t in ('aging','aging traits','aging (time to event)') and m == {'aging'}:
        return 'Combined aging traits','General or multivariate aging phenotype'
    return '', 'Outside scope: disease, cosmetic/organ-specific trait, cognition, component biomarker or unrelated phenotype'

candidate = json.loads((SRC / 'gwas_aging_search.json').read_text())
search_pattern=re.compile(r'\b(longevity|lifespan|life span|healthspan|health span|aging|ageing|epigenetic age|biological age|age at death)\b',re.I)
candidate=[r for r in candidate if search_pattern.search(r['DISEASE/TRAIT']+' '+r['MAPPED_TRAIT'])]
trait_audit = Counter(); gwas=[]; unresolved_gwas=[]
for r in candidate:
    category, reason=trait_category(r)
    trait_audit[(r['DISEASE/TRAIT'],r['MAPPED_TRAIT'],category or 'Excluded',reason)]+=1
    if not category: continue
    p=float(r['P-VALUE']) if r['P-VALUE'] else None
    included=p is not None and p <= 5e-8
    assoc_id=f'GW{len(gwas)+1:04d}'
    gene_modes=defaultdict(set)
    for field in ('REPORTED GENE(S)','MAPPED_GENE'):
        # Expand explicit gene-pair shorthand; preserve the unedited value in the association table.
        value=r[field].replace('CHRNA3/5','CHRNA3, CHRNA5').replace('NOL4L/-DT','NOL4L, NOL4L-DT')
        for token in re.split(r'\s+-\s+|\s*,\s*|\s*;\s*|\s*/\s*',value):
            token=token.strip()
            if not token: continue
            gene=resolve(token,source='GWAS '+field)
            if gene:
                gene_modes[gene].add('Author reported' if field=='REPORTED GENE(S)' else 'Catalog flanking' if ' - ' in r[field] else 'Catalog overlapping')
            elif included and token.upper() not in ('NR','INTERGENIC','NONE','NA','N/A','-'):
                unresolved_gwas.append({'Association':assoc_id,'Gene token':token,'Field':field,'Study':r['STUDY ACCESSION'],'Reason':'No unique HGNC match'})
    # ID mapping recovers Ensembl-only or obsolete displayed names without guessing aliases.
    for field in ('SNP_GENE_IDS','UPSTREAM_GENE_ID','DOWNSTREAM_GENE_ID'):
        for token in re.split(r'[;,\s]+',r[field]):
            if token and token not in ('NR','NA'):
                gene=resolve(token if token.startswith('ENSG') else '',token if token.isdigit() else '',source='GWAS gene ID')
                if gene: gene_modes[gene].add('Catalog overlapping' if field=='SNP_GENE_IDS' else 'Catalog flanking')
    url='https://www.ebi.ac.uk/gwas/studies/'+r['STUDY ACCESSION']
    gwas.append({'Association':assoc_id,'Trait category':category,'Reported trait':r['DISEASE/TRAIT'],
        'Included at P <= 5e-8':'Yes' if included else 'No','P value':p,'Human genes':joined(gene_modes),
        'Reported genes (original)':r['REPORTED GENE(S)'],'Mapped genes (original)':r['MAPPED_GENE'],
        'Variant':r['SNPS'],'Risk or effect allele':r['STRONGEST SNP-RISK ALLELE'],
        'Effect (source)':r['OR or BETA'],'CI / effect units':r['95% CI (TEXT)'],
        'P value context':r['P-VALUE (TEXT)'],'Study accession':r['STUDY ACCESSION'],
        'PMID':r['PUBMEDID'],'Paper':r['STUDY'],'Publication date':r['DATE'],
        'Initial sample':r['INITIAL SAMPLE SIZE'],'Replication sample':r['REPLICATION SAMPLE SIZE'],
        'Source URL':url})
    if included:
        for gene,modes in gene_modes.items():
            add(gene,'GWAS',category,assoc_id,'Homo sapiens',gene,r['DISEASE/TRAIT'],
                'Gene expression / perturbation direction not inferred',joined(modes),
                'Variant association; gene nominated by reported/proximity mapping, not established causality. '+r['P-VALUE (TEXT)'],
                url,pmid=r['PUBMEDID'],p=p,study=r['STUDY ACCESSION'],variant=r['SNPS'])

# Gene-level matrix. These source flags overlap and must not be added as independent evidence.
gene_evidence=defaultdict(list)
for r in evidence: gene_evidence[r['Human gene']].append(r)
master=[]
for gene,rows in sorted(gene_evidence.items()):
    h=by_symbol[gene]; source_set={r['Resource'] for r in rows}; types={r['Evidence type'] for r in rows}
    has=lambda t:'Yes' if t in types else 'No'
    gp=[r['P value'] for r in rows if r['Resource']=='GWAS']
    master.append({'Human gene':gene,'Data resources':'; '.join(s for s in ('Experiment','GenAge','CellAge','GWAS') if s in source_set),
        **{s:'Yes' if s in source_set else 'No' for s in ('Experiment','GenAge','CellAge','GWAS')},
        'Resource count':len(source_set),'Gene name':h['name'],'HGNC ID':h['hgnc_id'],'Entrez ID':h['entrez_id'],
        'Gene type':h['locus_type'],'GenAge human':has('Human candidate'),
        'GenAge animal homolog':has('Animal model homolog'),'GenAge fungal homolog':has('Fungal model homolog'),
        'CellAge manipulation':has('Senescence manipulation'),'CellAge expression':has('Senescence expression'),
        'GWAS longevity':has('Longevity / lifespan'),'GWAS healthspan':has('Healthspan / healthy aging'),
        'GWAS biological age':has('Biological / epigenetic aging'),'GWAS combined traits':has('Combined aging traits'),
        'Minimum GWAS P':min(gp) if gp else None,
        'Experiment organisms':joined(r['Tested organism'] for r in rows if r['Resource']=='Experiment'),
        'CellAge gene effect':joined(r['Direction in source'] for r in rows if r['Evidence type']=='Senescence manipulation'),
        'GenAge model influence':joined(r['Direction in source'] for r in rows if r['Resource']=='GenAge' and 'homolog' in r['Evidence type']),
        'Evidence rows':len(rows)})

traits=[{'Reported trait':t,'Mapped trait':m,'Category':c,'Selection reason':why,'Association rows':n}
        for (t,m,c,why),n in sorted(trait_audit.items())]
for records in (normalization,unresolved,mapping_audit):
    seen=set(); records[:]=[r for r in records if not (tuple(r.values()) in seen or seen.add(tuple(r.values())))]
stats={
    'human_genes':len(master),'source_counts':{s:sum(r[s]=='Yes' for r in master) for s in ('Experiment','GenAge','CellAge','GWAS')},
    'subtype_counts':{t:len({r['Human gene'] for r in evidence if r['Evidence type']==t}) for t in sorted({r['Evidence type'] for r in evidence})},
    'at_least_two_resources':sum(r['Resource count']>=2 for r in master),
    'evidence_rows':len(evidence),'experiment_records':len(experiments),
    'experiment_gene_species_pairs':len({(r['Scientific name'],r['Gene in tested organism']) for r in experiments}),
    'experiment_unmapped_records':[r['Record ID'] for r in experiments if not r['Human candidates']],
    'genage_human_source_records':len(human_genage),'genage_model_source_records':len(models),
    'genage_model_mapped':sum(bool(r['Human candidates']) for r in models),
    'cellage_manipulation_records':len(cellage),'cellage_source_genes':len({r['Entrez ID'] for r in cellage}),
    'cellage_expression_records':len(signatures),
    'gwas_in_scope_associations':len(gwas),'gwas_significant_associations':sum(r['Included at P <= 5e-8']=='Yes' for r in gwas),
    'gwas_significant_studies':len({r['Study accession'] for r in gwas if r['Included at P <= 5e-8']=='Yes'}),
    'gwas_significant_publications':len({r['PMID'] for r in gwas if r['Included at P <= 5e-8']=='Yes'}),
    'gwas_unmapped_significant_associations':sum(r['Included at P <= 5e-8']=='Yes' and not r['Human genes'] for r in gwas),
    'gwas_unresolved_tokens':len(unresolved_gwas),
    'gwas_release':json.loads((SRC/'gwas_metadata.json').read_text())['data_release_date'],
}
data=dict(master=master,evidence=evidence,models=models,experiments=experiments,gwas=gwas,traits=traits,stats=stats)
(ROOT/'reference.json').write_text(json.dumps(data,indent=2))
(ROOT/'summary.json').write_text(json.dumps(stats,indent=2))
for name,rows in [('aging_longevity_genes.tsv',master),('evidence_details.tsv',evidence),('genage_model_records.tsv',models),
                  ('experimental_records.tsv',experiments),('gwas_associations.tsv',gwas),('gwas_trait_selection.tsv',traits),
                  ('ortholog_mapping_audit.tsv',mapping_audit),('symbol_normalization.tsv',normalization),
                  ('unresolved_identifiers.tsv',unresolved),('unresolved_gwas_tokens.tsv',unresolved_gwas)]:
    write_tsv(name,rows)
(ROOT/'all_human_gene_symbols.txt').write_text('\n'.join(r['Human gene'] for r in master)+'\n')
for source in ('Experiment','GenAge','CellAge','GWAS'):
    (ROOT/f'{source.lower()}_human_genes.txt').write_text('\n'.join(r['Human gene'] for r in master if r[source]=='Yes')+'\n')
for field,fn in [('CellAge manipulation','cellage_manipulation_genes.txt'),('CellAge expression','cellage_expression_genes.txt'),
                 ('GWAS longevity','gwas_longevity_genes.txt'),('GWAS healthspan','gwas_healthspan_genes.txt'),
                 ('GWAS biological age','gwas_biological_age_genes.txt'),('GWAS combined traits','gwas_combined_traits_genes.txt')]:
    (ROOT/fn).write_text('\n'.join(r['Human gene'] for r in master if r[field]=='Yes')+'\n')
manifest={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in SRC.iterdir() if p.is_file()}
(ROOT/'source_checksums.json').write_text(json.dumps(manifest,indent=2))
assert len(master)==len({r['Human gene'] for r in master})
assert all(r['Human gene'] in by_symbol for r in master)
assert len(models)==2205 and len(human_genage)==307 and len({r['Entrez ID'] for r in cellage})==866
assert all(r['P value']<=5e-8 for r in evidence if r['Resource']=='GWAS')
assert all(not re.search(r'imaging|reproductive lifespan',r['Phenotype'],re.I) for r in evidence if r['Resource']=='GWAS')
print(json.dumps(stats,indent=2))
