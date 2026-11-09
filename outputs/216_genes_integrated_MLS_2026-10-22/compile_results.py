"""Integrate the 216 regulatory hits with existing aging evidence and assay fits."""
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
REF = ROOT / 'outputs/aging_longevity_reference_2026-10-22'
REG = ROOT / 'outputs/regulatory_windows_2026-10-22'
ACT = ROOT / 'outputs/activity_omics_comparison_2026-10-20'
WEB = ROOT / 'cross_species_visual_tool/web_app'
sys.path[:0] = [str(ROOT / 'output/activity_omics_comparison_2026-10-20'), str(WEB), str(WEB / '.deps')]
import numpy as np
from scipy.stats import t, false_discovery_control
import quantify_activity as qa
from data_paths import workbook_path


def read(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter='\t'))


def write(name, records):
    with (OUT / name).open('w') as h:
        w = csv.DictWriter(h, fieldnames=list(records[0]), delimiter='\t')
        w.writeheader()
        w.writerows(records)


sig = read(REG / 'significant_genes_all_windows.tsv')
by_gene = defaultdict(list)
for r in sig:
    assert float(r['BH_FDR_window']) < .05
    by_gene[r['gene']].append(r)
assert len(sig) == 300 and len(by_gene) == 216
assert all(len({r['direction'] for r in rows}) == 1 for rows in by_gene.values())
genes = sorted(by_gene)

hgnc = read(REF / 'sources/hgnc_complete_set.tsv')
approved = {r['symbol'].upper(): r['symbol'] for r in hgnc}
aliases = defaultdict(set)
for r in hgnc:
    for a in (r['alias_symbol'] + '|' + r['prev_symbol']).split('|'):
        if a:
            aliases[a.upper()].add(r['symbol'])


def canonical(g):
    key = g.upper()
    return approved.get(key) or (next(iter(aliases[key])) if len(aliases[key]) == 1 else key)


assert len({canonical(g) for g in genes}) == 216
reference = {r['Human gene']: r for r in read(REF / 'aging_longevity_genes.tsv')}
labels = {'tss2kb': 'Gene activity: TSS ±2 kb', 'body': 'Gene activity: body + 2 kb',
          'rna': 'Transcriptome', 'protein': 'Proteome'}
fits, canonical_fit = {}, {}
for m in labels:
    rows = read(ACT / f'mass_adjusted/all_available/{m}_correlations.tsv')
    assert np.allclose(false_discovery_control([float(r['p']) for r in rows]),
                       [float(r['q']) for r in rows], rtol=1e-12, atol=1e-14)
    fits[m] = {r['gene']: r for r in rows}
    canonical_fit[m] = defaultdict(list)
    for r in rows:
        canonical_fit[m][canonical(r['gene'])].append(r['gene'])


def lookup(g, m):
    if g in fits[m]:
        return g, 'Exact source symbol'
    c = canonical(g)
    if c in fits[m]:
        return c, 'Approved human symbol'
    matches = canonical_fit[m][c]
    if len(matches) == 1:
        return matches[0], 'Unique HGNC previous/alias symbol'
    return None, 'Ambiguous symbol match' if len(matches) > 1 else 'No eligible fitted result'


# Check source versions before joining the already computed genome-wide fits.
manifest = json.loads((REG / 'manifest.json').read_text())
assert hashlib.sha256((WEB / 'data/species_mls.json').read_bytes()).hexdigest() == manifest['traits_sha256']
assert (WEB / 'data/gene_activity.sqlite').stat().st_mtime_ns == manifest['database_mtime_ns']
for r in manifest['source_files']:
    p = Path(r['path'])
    assert p.stat().st_size == r['size'] and p.stat().st_mtime_ns == r['mtime_ns'], str(p)

specs = {s[0]: s for s in qa.SELECTED}
old_traits = {r['species']: r for r in read(ACT / 'species_inclusion.tsv')}
assert set(old_traits) == set(specs)
cache_validations = []
selected_cache = {}
for sid, s in specs.items():
    assert float(old_traits[sid]['MLS_years']) == qa.LIFE[s[2]]['years']
    assert float(old_traits[sid]['body_mass_g']) == qa.LIFE[s[2]]['body_mass_g']
    loci = qa.read_loci(sid)
    files = sorted(p for p in (ROOT / 'cross_species_visual_tool' / s[4]).glob('*.bw')
                   if p.name.lower().startswith(s[5].lower()) and
                   p.stem.removesuffix('.cpm').removesuffix('_scaled') not in qa.EXCLUDED_SAMPLES)
    fp_input = dict(version=1, loci_sha256=hashlib.sha256(json.dumps(loci).encode()).hexdigest(),
                    methods=qa.METHODS, files=[dict(path=str(p.resolve()), size=p.stat().st_size,
                    mtime_ns=p.stat().st_mtime_ns) for p in files])
    fp = hashlib.sha256(json.dumps(fp_input, sort_keys=True).encode()).hexdigest()
    assert fp == json.loads((qa.CACHE / f'{sid}.json').read_text())['fingerprint'], sid
    cache_validations.append(dict(species=sid, fingerprint=fp))
    with gzip.open(qa.CACHE / f'{sid}.tsv.gz', 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['gene'] in by_gene and r['method'] in ('tss2kb', 'body'):
                if r['region_eligible'] == '1' and float(r['median_activity']) > 0:
                    selected_cache[r['method'], r['gene'], sid] = r
print('Verified regulatory inputs, activity caches and species traits.', flush=True)

source_settings = json.loads((ACT / 'settings.json').read_text())
with sqlite3.connect(f'file:{WEB}/data/liver_omics.sqlite?mode=ro', uri=True) as con:
    omics_sources = json.loads(con.execute("SELECT value FROM metadata WHERE key='sources'").fetchone()[0])
assert omics_sources == source_settings['omics_sources']
for r in omics_sources.values():
    p = workbook_path(r['file'])
    assert p.stat().st_size == r['size'] and p.stat().st_mtime_ns == r['mtime_ns'], r['file']

required = set(genes)
for g in genes:
    for m in labels:
        symbol, _ = lookup(g, m)
        if symbol:
            required.add(symbol)
inputs = defaultdict(dict)
lengths = {}
with gzip.open(ACT / 'analysis_input_species_medians.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        if r['method'] in labels and r['gene'] in required:
            m, g, s = r['method'], r['gene'], r['species']
            inputs[m, g][s] = float(r['analysis_value'])
            if m == 'body':
                lengths[g, s] = float(r['log10_gene_length_kb'])
for (m, g, s), r in selected_cache.items():
    assert math.isclose(inputs[m, g][s], math.log10(float(r['median_activity'])), abs_tol=1e-12)

# Verify selected RNA/protein medians against the current workbook-derived index.
latin_to_sid = {s[2]: sid for sid, s in specs.items()}
latin_to_sid['Ellobius lutescens'] = 'Bramus_lutescens'
raw = defaultdict(list)
with sqlite3.connect(f'file:{WEB}/data/liver_omics.sqlite?mode=ro', uri=True) as con:
    placeholders = ','.join('?' for _ in required)
    query = f'SELECT assay,symbol,species,log,normalized,norm_on FROM measurements WHERE symbol IN ({placeholders}) AND NOT (assay="rna" AND sample="BMR_Turk_Lu")'
    for assay, g, latin, log, normalized, norm_on in con.execute(query, sorted(required)):
        if latin not in latin_to_sid:
            continue
        value = log if assay == 'rna' else normalized
        if value is None or not math.isfinite(value):
            continue
        if assay == 'protein':
            assert norm_on == 'log.quant'
        raw[assay, g, latin_to_sid[latin]].append(value)
for (m, g), rows in inputs.items():
    if m in ('rna', 'protein'):
        for sid, value in rows.items():
            assert math.isclose(float(np.median(raw[m, g, sid])), value, rel_tol=1e-12, abs_tol=1e-12), (m, g, sid)
print('Verified selected RNA/protein medians against the current index.', flush=True)

details, overview, audit = [], [], []
for g in genes:
    c = canonical(g)
    ref = reference.get(c, {})
    windows = sorted(int(r['window_kb']) for r in by_gene[g])
    representative = min(by_gene[g], key=lambda r: (float(r['BH_FDR_window']), float(r['p_value']), int(r['window_kb'])))
    direction = representative['direction']
    row = dict(Gene=g, Human_symbol=c, Significant_windows_kb='; '.join(map(str, windows)),
               Experiment=ref.get('Experiment', 'No'), GenAge=ref.get('GenAge', 'No'),
               CellAge=ref.get('CellAge', 'No'), GWAS=ref.get('GWAS', 'No'),
               Regulatory=f'{direction} **')
    details.append(dict(Gene=g, Human_symbol=c, Method='Regulatory peak count', Direction=direction,
        Partial_r=float(representative['partial_r']), P=float(representative['p_value']),
        BH_FDR=float(representative['BH_FDR_window']), N_species=int(representative['n_species']),
        MLS_beta=float(representative['MLS_coefficient']), DF=int(representative['df']),
        Status='FDR < 0.05', Representative_window_kb=int(representative['window_kb']),
        Matched_source_symbol=g, Symbol_matching='Exact source symbol',
        Species_IDs=representative['species'], FDR_family='All tested genes within this window'))
    for m, label in labels.items():
        key, match_rule = lookup(g, m)
        saved = fits[m].get(key)
        if saved:
            ids = saved['species'].split(',')
            assert set(ids) == set(inputs[m, key]), (m, key, 'species')
            x = np.array([math.log10(qa.LIFE[specs[s][2]]['years']) for s in ids])
            y = np.array([inputs[m, key][s] for s in ids])
            mass = np.array([math.log10(qa.LIFE[specs[s][2]]['body_mass_g']) for s in ids])
            columns = [np.ones(len(ids)), x, mass]
            if m == 'body':
                columns.append(np.array([lengths[key, s] for s in ids]))
            design = np.column_stack(columns)
            beta, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
            assert rank == design.shape[1]
            df = len(ids) - int(rank)
            se = math.sqrt(float(np.sum((y-design@beta)**2))/df * np.linalg.inv(design.T@design)[1, 1])
            stat = beta[1] / se
            p = float(2*t.sf(abs(stat), df))
            r = float(stat / math.sqrt(stat**2 + df))
            for name, val in [('r', r), ('p', p), ('MLS_beta', float(beta[1]))]:
                assert math.isclose(val, float(saved[name]), rel_tol=1e-7, abs_tol=1e-10), (m, g, name)
            q = float(saved['q'])
            direction = 'Positive' if r > 0 else 'Negative'
            mark = ' **' if q < .05 else ' *' if p < .05 else ''
            row[m] = direction + mark
            status = 'FDR < 0.05' if q < .05 else 'Nominal P < 0.05 only' if p < .05 else 'Not significant'
            details.append(dict(Gene=g, Human_symbol=c, Method=label, Direction=direction,
                Partial_r=r, P=p, BH_FDR=q, N_species=len(ids), MLS_beta=float(beta[1]), DF=df,
                Status=status, Representative_window_kb=None, Matched_source_symbol=key,
                Symbol_matching=match_rule, Species_IDs=','.join(ids),
                FDR_family=f'All {len(fits[m]):,} tested genes in this method'))
            audit.append(dict(gene=g, method=m, source_symbol=key, P_OLS=p, P_source=float(saved['p'])))
        else:
            ids = sorted(inputs[m, g])
            reason = 'No eligible species measurements' if not ids else f'Only {len(ids)} eligible species; minimum 5' if len(ids) < 5 else 'No valid model fit in source results'
            if match_rule.startswith('Ambiguous'):
                reason = match_rule
            row[m] = 'Not tested'
            details.append(dict(Gene=g, Human_symbol=c, Method=label, Direction='Not tested',
                Partial_r=None, P=None, BH_FDR=None, N_species=len(ids), MLS_beta=None, DF=None,
                Status=reason, Representative_window_kb=None, Matched_source_symbol=g,
                Symbol_matching=match_rule, Species_IDs=','.join(ids),
                FDR_family=f'All {len(fits[m]):,} tested genes in this method'))
    row['Aging_resources'] = ref.get('Data resources', 'Not found in compiled reference')
    overview.append(row)

assert len(overview) == 216 and len(details) == 1080
assert sum(r['Aging_resources'] != 'Not found in compiled reference' for r in overview) == 46
evidence = [r for r in read(REF / 'evidence_details.tsv') if r['Human gene'] in {canonical(g) for g in genes}]
summary = []
for m in ['Regulatory peak count'] + list(labels.values()):
    rows = [r for r in details if r['Method'] == m]
    valid = [r for r in rows if r['P'] is not None]
    summary.append(dict(Method=m, Positive=sum(r['Direction']=='Positive' for r in rows),
        Negative=sum(r['Direction']=='Negative' for r in rows), Not_tested=216-len(valid),
        Nominal_positive=sum(r['P'] < .05 and r['Direction']=='Positive' for r in valid),
        Nominal_negative=sum(r['P'] < .05 and r['Direction']=='Negative' for r in valid),
        FDR_positive=sum(r['BH_FDR'] < .05 and r['Direction']=='Positive' for r in valid),
        FDR_negative=sum(r['BH_FDR'] < .05 and r['Direction']=='Negative' for r in valid),
        Min_species=min(r['N_species'] for r in valid), Max_species=max(r['N_species'] for r in valid)))
sources = Counter(s for r in overview for s in ['Experiment','GenAge','CellAge','GWAS'] if r[s]=='Yes')
notes = [
    {'Item':'Scope','Details':'216 unique genes significant in at least one regulatory TSS window: ±20, ±50, ±100, ±200, ±500 or ±1000 kb; 300 significant gene-window pairs. Each window is a flank: total interval width is twice the stated value.'},
    {'Item':'Direction symbols','Details':'Positive/Negative denotes the sign of the MLS partial correlation. ** means BH FDR < 0.05; * means nominal P < 0.05 only. Unmarked directions are not nominally significant. Not tested means no eligible fitted result; it is not a negative or zero correlation.'},
    {'Item':'One regulatory direction','Details':'All significant windows have the same direction within each gene: 215 genes positive, FBXO25 negative. Overview lists every significant window and shows the direction once. For the detailed effect/P/FDR/N, select the window with smallest within-window FDR; ties use P then the smaller window. This representative selection is descriptive, not an additional significance test.'},
    {'Item':'Regulatory measure','Details':'log10(species median of sample-specific ATAC peak counts overlapping the TSS window), modeled against log10(MLS) with log10(body mass) as covariate. No length adjustment.'},
    {'Item':'Gene activity measures','Details':'Unweighted integrated bigWig CPM·bp across the complete region, not restricted to peaks: TSS ±2 kb (4 kb total) or whole gene body + 2 kb upstream promoter. Use log10 of the species median. Both adjust for log10(body mass); gene body also adjusts for log10(gene-body length).'},
    {'Item':'Transcriptome','Details':'Liver species median of source log.tpm; no additional logarithm. Correlation with log10(MLS) adjusts for log10(body mass). BMR_Turk_Lu excluded. Mouse liver supplement retained. Source RNA annotation labels ASM as Acomys cahirinus but common field says Acomys russatus; current source mapping is retained.'},
    {'Item':'Proteome','Details':'Liver species median of source quant.norm normalized on log.quant; no additional logarithm. Correlation with log10(MLS) adjusts for log10(body mass).'},
    {'Item':'Species coverage','Details':'Each method uses all eligible species for each gene; species are NOT forced to match across assays. Human, rhesus macaque and bowhead whale excluded. Missing/lower-bound MLS or missing mass excluded. Minimum 5 species per fit. Correlation details retain N and exact species IDs; sparse proteome coverage remains Not tested.'},
    {'Item':'ATAC region eligibility','Details':'Choose the locus with highest species-median activity separately for each method before filtering, without fallback. Exclude nonpositive values for log scales and regions with <50% available genomic span. Availability uses contig boundaries, not ATAC coverage or internal assembly gaps.'},
    {'Item':'FDR families','Details':f'BH correction is retained from all tested genes within each regulatory window and from the full genome-wide assay analysis: TSS ±2 kb {len(fits["tss2kb"]):,}; gene body {len(fits["body"]):,}; RNA {len(fits["rna"]):,}; protein {len(fits["protein"]):,}. It is not recalculated over these 216 selected genes. Joint-across-window FDR is not the primary selection criterion.'},
    {'Item':'Aging evidence flags','Details':'Yes means present in the previously compiled source-labelled reference; No means not found in this collection, not proof that no aging study exists. 46 genes have at least one source; 170 have no match. GenAge includes human candidates and model-organism homologs; CellAge includes manipulation and expression evidence; GWAS assignments include reported/proximity genes and are not necessarily causal.'},
    {'Item':'Experimental evidence','Details':'Experiment flags the previously curated animal lifespan-extension perturbation collection, not all experimentally studied genes in GenAge or CellAge. GenAge/CellAge source names are preserved. Cross-species ATAC directions are not inferred perturbation effects; source-organism evidence is not proof of an effect in humans.'},
    {'Item':'Symbol handling','Details':'Original regulatory symbols are retained with HGNC-approved human symbols alongside. MFSD10 → SLC75A1; B3GNTL1 → QTGAL; TNFAIP8L1 → TIPE1. DENND10 RNA/protein results are recovered via its unique previous name FAM45A. Exact assay symbol is preferred, then approved symbol, then a unique HGNC alias; ambiguous matches are not forced.'},
    {'Item':'Interpretation','Details':'These are observational cross-species associations without phylogenetic, sequencing-depth or batch correction. Species sets and selected loci can differ. Gene-symbol matching does not verify one-to-one orthology. Follow-up assays use overlapping samples, so concordance is exploratory, not independent replication.'},
    {'Item':'Verification','Details':f'Checked current regulatory source metadata, activity fingerprints, trait values and RNA/protein workbook provenance. Selected RNA/protein medians were checked against the current index. All {len(audit)} available gene-activity/RNA/protein fits were independently checked using full-model OLS; genome-wide BH values were verified.'},
    {'Item':'Regulatory results','Details':str(REG / 'significant_genes_all_windows.tsv')},
    {'Item':'Activity / expression results','Details':str(ACT / 'mass_adjusted/all_available')},
    {'Item':'Aging reference and primary links','Details':str(REF / 'aging_longevity_genes.tsv') + '; evidence_details.tsv contains the study/database URLs reproduced on Aging evidence.'},
    {'Item':'Species traits','Details':str(WEB / 'data/species_mls.json')},
]
write('216_genes_overview.tsv', overview)
write('correlation_details.tsv', details)
write('aging_evidence.tsv', evidence)
write('direction_summary.tsv', summary)
write('methods.tsv', notes)
data = dict(overview=overview, details=details, evidence=evidence, summary=summary, notes=notes,
            source_counts=dict(sources), supported_genes=46, unsupported_genes=170)
(OUT / 'results.json').write_text(json.dumps(data, indent=2))
input_files = [REG / 'significant_genes_all_windows.tsv', REF / 'aging_longevity_genes.tsv',
               REF / 'evidence_details.tsv', ACT / 'settings.json'] + [
               ACT / f'mass_adjusted/all_available/{m}_correlations.tsv' for m in labels]
(OUT / 'validation.json').write_text(json.dumps(dict(activity_cache=cache_validations, OLS=audit,
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in input_files}), indent=2))
print(json.dumps(dict(summary=summary, evidence_rows=len(evidence), source_counts=dict(sources),
    recovered_symbols=[r for r in details if r['Symbol_matching']=='Unique HGNC previous/alias symbol']), indent=2))
