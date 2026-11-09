"""Join corrected/uncorrected regulatory hits with the existing assay/evidence data."""
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, t

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
PREV = ROOT / 'outputs/216_genes_integrated_MLS_2026-10-22'
REF = ROOT / 'outputs/aging_longevity_reference_2026-10-22'
REG = ROOT / 'outputs/regulatory_windows_2026-10-22'
RAW = ROOT / 'outputs/regulatory_windows_no_mass_2026-10-23'
ACT = ROOT / 'outputs/activity_omics_comparison_2026-10-20'
WINDOWS = [20, 50, 100, 200, 500, 1000]
LABELS = {'tss2kb': 'Gene activity TSS ±2 kb', 'body': 'Gene activity body + 2 kb',
          'rna': 'Liver transcriptome', 'protein': 'Liver proteome'}


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as h:
        return list(csv.DictReader(h, delimiter='\t'))


def write(name, rows):
    with (OUT/name).open('w', newline='') as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


previous = json.loads((PREV/'results.json').read_text())
old_audit = json.loads((PREV/'validation.json').read_text())
for path, checksum in old_audit['source_sha256'].items():
    assert sha(Path(path)) == checksum, f'Changed previously validated input: {path}'
corrected = read(REG/'significant_genes_all_windows.tsv')
uncorrected = read(RAW/'significant_genes_all_windows.tsv')
sig = {'Corrected': defaultdict(list), 'Uncorrected': defaultdict(list)}
for model, records in [('Corrected', corrected), ('Uncorrected', uncorrected)]:
    for r in records:
        assert float(r['BH_FDR_window']) < .05
        sig[model][r['gene']].append(r)
gene_sets = {m: set(records) for m, records in sig.items()}
genes = sorted(set.union(*gene_sets.values()))
assert len(genes) == 330
assert len(gene_sets['Corrected'] & gene_sets['Uncorrected']) == 143

hgnc = read(REF/'sources/hgnc_complete_set.tsv')
approved = {r['symbol'].upper(): r['symbol'] for r in hgnc}
aliases = defaultdict(set)
for r in hgnc:
    for a in (r['alias_symbol']+'|'+r['prev_symbol']).split('|'):
        if a:
            aliases[a.upper()].add(r['symbol'])


def canonical(g):
    return approved.get(g.upper()) or (next(iter(aliases[g.upper()])) if len(aliases[g.upper()]) == 1 else g.upper())


assert len({canonical(g) for g in genes}) == 330
reference = {r['Human gene']: r for r in read(REF/'aging_longevity_genes.tsv')}
fits, fit_canonical, input_paths = {}, {}, []
for method in LABELS:
    path = ACT/f'mass_adjusted/all_available/{method}_correlations.tsv'
    input_paths.append(path)
    rows = read(path)
    assert np.allclose(false_discovery_control([float(r['p']) for r in rows]),
                       [float(r['q']) for r in rows], rtol=1e-12, atol=1e-14)
    fits[method] = {r['gene']: r for r in rows}
    fit_canonical[method] = defaultdict(list)
    for r in rows:
        fit_canonical[method][canonical(r['gene'])].append(r['gene'])


def lookup(g, method):
    if g in fits[method]:
        return g, 'Exact source symbol'
    c = canonical(g)
    if c in fits[method]:
        return c, 'Approved human symbol'
    matches = fit_canonical[method][c]
    return (matches[0], 'Unique HGNC previous/alias symbol') if len(matches) == 1 else (None, 'Ambiguous symbol match' if matches else 'No eligible fitted result')


required = set(genes)
for g in genes:
    for m in LABELS:
        key, _ = lookup(g, m)
        if key:
            required.add(key)
inputs, lengths = defaultdict(dict), {}
with gzip.open(ACT/'analysis_input_species_medians.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        if r['method'] in LABELS and r['gene'] in required:
            m, g, sid = r['method'], r['gene'], r['species']
            inputs[m, g][sid] = float(r['analysis_value'])
            if m == 'body':
                lengths[g, sid] = float(r['log10_gene_length_kb'])
traits = {r['species']: r for r in read(ACT/'species_inclusion.tsv')}
current_life = json.loads((ROOT/'cross_species_visual_tool/web_app/data/species_mls.json').read_text())
for r in read(REG/'species.tsv'):
    if r['species'] in traits:
        sid = r['species']
        assert float(traits[sid]['MLS_years']) == current_life[r['scientific_name']]['years']
        assert float(traits[sid]['body_mass_g']) == current_life[r['scientific_name']]['body_mass_g']

overview, details, audit = [], [], []
for g in genes:
    c = canonical(g)
    ref = reference.get(c, {})
    membership = 'Both' if all(g in s for s in gene_sets.values()) else 'Corrected only' if g in gene_sets['Corrected'] else 'Uncorrected only'
    row = dict(Gene=g, Human_symbol=c, Regulatory_significance=membership)
    for model in ['Corrected', 'Uncorrected']:
        records = sig[model].get(g, [])
        row[f'{model}_windows_kb'] = '; '.join(map(str, sorted(int(r['window_kb']) for r in records))) if records else 'None'
    for source in ['Experiment', 'GenAge', 'CellAge', 'GWAS']:
        row[source] = ref.get(source, 'No')
    row['Aging_resources'] = ref.get('Data resources', 'Not found in compiled reference')
    for model in ['Corrected', 'Uncorrected']:
        directions = {r['direction'] for r in sig[model].get(g, [])}
        row[f'Regulatory_{model.lower()}'] = next(iter(directions))+' **' if len(directions) == 1 else 'Mixed **' if directions else 'Not FDR-significant'
    for m, label in LABELS.items():
        key, match = lookup(g, m)
        saved = fits[m].get(key)
        ids = saved['species'].split(',') if saved else sorted(inputs[m, g])
        if saved:
            assert set(ids) == set(inputs[m, key])
            x = np.log10([float(traits[s]['MLS_years']) for s in ids])
            y = np.array([inputs[m, key][s] for s in ids])
            mass = np.log10([float(traits[s]['body_mass_g']) for s in ids])
            columns = [np.ones(len(ids)), x, mass]
            if m == 'body':
                columns.append(np.array([lengths[key, s] for s in ids]))
            X = np.column_stack(columns)
            beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
            assert rank == X.shape[1]
            df = len(ids)-rank
            se = math.sqrt(float(np.sum((y-X@beta)**2))/df*np.linalg.inv(X.T@X)[1, 1])
            stat = beta[1]/se
            p = float(2*t.sf(abs(stat), df))
            r = float(stat/math.sqrt(stat*stat+df))
            for name, value in [('p', p), ('r', r), ('MLS_beta', float(beta[1]))]:
                assert math.isclose(value, float(saved[name]), rel_tol=1e-7, abs_tol=1e-10), (g, m, name)
            q = float(saved['q'])
            direction = 'Positive' if r > 0 else 'Negative' if r < 0 else 'Zero'
            row[m] = direction+(' **' if q < .05 else ' *' if p < .05 else '')
            status = 'FDR < 0.05' if q < .05 else 'Nominal P < 0.05 only' if p < .05 else 'Not significant'
            audit.append(dict(gene=g, method=m, source_symbol=key, p_difference=abs(p-float(saved['p']))))
        else:
            r = p = q = df = None
            direction = row[m] = 'Not tested'
            status = match if match.startswith('Ambiguous') else 'No eligible species measurements' if not ids else f'Only {len(ids)} eligible species; minimum 5' if len(ids) < 5 else 'No valid source fit'
        details.append(dict(Gene=g, Human_symbol=c, Method=label, Body_mass_correction='Yes',
            Covariates='log10(body mass), log10(gene-body length)' if m == 'body' else 'log10(body mass)',
            Direction=direction, Partial_r=float(saved['r']) if saved else None, P=float(saved['p']) if saved else None,
            BH_FDR=q, N_species=len(ids), DF=int(df) if df is not None else None,
            MLS_coefficient=float(saved['MLS_beta']) if saved else None, Status=status,
            Source_gene_symbol=key or g, Symbol_matching=match,
            FDR_family=f'All {len(fits[m]):,} tested genes in this method', Species_IDs=','.join(ids)))
    overview.append(row)
print(f'Joined {len(overview)} genes and independently verified {len(audit)} assay fits.', flush=True)

evidence = [r for r in read(REF/'evidence_details.tsv') if r['Human gene'] in {canonical(g) for g in genes}]
by_evidence = defaultdict(list)
for e in evidence:
    resource, typ, direction = e['Resource'], e['Evidence type'], e['Direction in source']
    role, scope = 'Unclassified', 'No directional perturbation evidence in this record'
    if resource == 'Experiment' and typ == 'Animal lifespan extension':
        if direction in ['Increase', 'Decrease']:
            role = 'Pro-longevity' if direction == 'Increase' else 'Anti-longevity'
        elif e['Source record'] == 'EXP040' and e['Original gene'] == 'Thor':
            role = 'Pro-longevity'
        scope = 'Lifespan extension intervention in '+e['Tested organism']
        if e['Source record'] == 'EXP040':
            scope += '; activated fly Thor transgene, not ordinary wild-type overexpression; human ortholog inference only'
    elif resource == 'GenAge' and direction in ['Pro-Longevity', 'Anti-Longevity']:
        role = 'Pro-longevity' if direction == 'Pro-Longevity' else 'Anti-longevity'
        scope = 'Source-organism label: '+e['Tested organism']+'; human homolog inference only'
    elif resource == 'CellAge' and typ == 'Senescence manipulation':
        if direction in ['Induces', 'Inhibits']:
            role = 'Pro-senescence' if direction == 'Induces' else 'Anti-senescence'
        scope = 'Human cellular senescence; not a demonstrated whole-animal aging/lifespan role'
    elif resource == 'GWAS':
        scope = 'Variant association; gene functional direction not established'
    elif typ == 'Senescence expression':
        scope = 'Expression association; pro-/anti-aging causality not established'
    elif direction == 'Necessary for fitness':
        scope = 'Fitness evidence; not reclassified as pro-longevity'
    e['Interpreted role'], e['Role scope'] = role, scope
    by_evidence[e['Human gene']].append(e)

roles = []
for row in overview:
    records = by_evidence[row['Human_symbol']]
    life = {r['Interpreted role'] for r in records if r['Interpreted role'] in ['Pro-longevity', 'Anti-longevity']}
    cell = {r['Interpreted role'] for r in records if r['Interpreted role'] in ['Pro-senescence', 'Anti-senescence']}
    life_role = 'Mixed / context-dependent' if len(life) > 1 else next(iter(life)) if life else 'Unknown'
    cell_role = 'Mixed / context-dependent' if len(cell) > 1 else next(iter(cell)) if cell else 'Unknown'
    pieces = []
    if len(life) > 1:
        pieces.append('Mixed longevity roles (model homologs)')
    elif life:
        experimental = [r for r in records if r['Resource']=='Experiment' and r['Interpreted role'] in life]
        if experimental:
            animals = {r['Tested organism'] for r in experimental}
            scope = 'mouse experiment + model evidence' if row['Gene']=='SIRT6' else 'mouse experiments' if animals=={'Mus musculus'} else 'fly activation experiment; homolog evidence' if animals=={'Drosophila melanogaster'} else 'animal experiment; homolog evidence'
        else:
            scope = 'model homolog evidence'
        pieces.append(next(iter(life))+f' ({scope})')
    if len(cell) > 1:
        pieces.append('Mixed senescence roles (cells)')
    elif cell:
        pieces.append(next(iter(cell))+' (cellular aging only)')
    role = '\n'.join(pieces or ['Unknown / not established' if records else 'Unknown (not in compiled reference)'])
    row['Reported_aging_longevity_role'] = role
    directional = [e for e in records if e['Interpreted role'] != 'Unclassified']
    basis = '; '.join(sorted({f"{e['Resource']}: {e['Original gene']} ({e['Tested organism']}), {e['Direction in source']}" for e in directional}))
    if not basis and records:
        basis = 'Existing association, human-candidate or fitness evidence does not establish a pro-/anti-aging role.'
    roles.append(dict(Gene=row['Gene'], Human_symbol=row['Human_symbol'], Reported_role=role,
        Longevity_role=life_role, Cellular_senescence_role=cell_role, Directional_evidence=bool(life or cell),
        Evidence_basis=basis, Interpretation='Roles describe the cited source context; homologs and cellular senescence are not established human lifespan effects.',
        Source_URLs=' ; '.join(sorted({e['Source URL'] for e in directional}))))

regulatory = []
for model, folder, r_field in [('Corrected', REG, 'partial_r'), ('Uncorrected', RAW, 'r')]:
    for kb in WINDOWS:
        path = folder/f'tss_{kb}kb/all_tested.tsv.gz'
        input_paths.append(path)
        source = read(path)
        assert np.allclose(false_discovery_control([float(r['p_value']) for r in source]),
                           [float(r['BH_FDR_window']) for r in source], rtol=1e-12, atol=1e-14)
        mapping = {r['gene']: r for r in source}
        for g in genes:
            r = mapping.get(g)
            regulatory.append(dict(Gene=g, Human_symbol=canonical(g), Model=model,
                Body_mass_correction='Yes' if model=='Corrected' else 'No', Window_flank_kb=kb,
                Full_width_kb=2*kb, Direction=('Positive' if float(r[r_field])>0 else 'Negative') if r else 'Not tested',
                Correlation_r=float(r[r_field]) if r else None, P=float(r['p_value']) if r else None,
                BH_FDR_window=float(r['BH_FDR_window']) if r else None,
                Significant_FDR05='Yes' if r and float(r['BH_FDR_window'])<.05 else 'No' if r else 'Not tested',
                N_species=int(r['n_species']) if r else None, DF=int(r['df']) if r else None,
                MLS_coefficient=float(r['MLS_coefficient']) if r else None,
                FDR_family=f'All {len(source):,} tested genes within this model and window',
                Species_IDs=r['species'] if r else ''))
regulatory.sort(key=lambda r:(r['Gene'],r['Window_flank_kb'],r['Model']))

# Preserve the previous 216-gene evidence, roles and all four assay directions exactly.
mapping = {r['Gene']:r for r in overview}
for old in previous['overview']:
    new = mapping[old['Gene']]
    for key in ['Human_symbol','Experiment','GenAge','CellAge','GWAS','Reported_aging_longevity_role','tss2kb','body','rna','protein']:
        assert new[key] == old[key], (old['Gene'],key,new[key],old[key])
    assert new['Corrected_windows_kb'] == old['Significant_windows_kb']
    assert new['Regulatory_corrected'] == old['Regulatory']
categories = dict(Counter(r['Regulatory_significance'] for r in overview))
assert categories == {'Both':143,'Corrected only':73,'Uncorrected only':114}
supported = sum(any(r[k]=='Yes' for k in ['Experiment','GenAge','CellAge','GWAS']) for r in overview)
notes = [
    dict(Item='Scope', Details='330 unique genes significant at BH FDR < 0.05 in any of six regulatory windows in either model. 216 corrected, 257 uncorrected: 143 shared, 73 corrected only and 114 uncorrected only.'),
    dict(Item='Correction status', Details='Both / Corrected only / Uncorrected only describes regulatory significance in at least one window. Corrected includes log10(body mass in g) as a covariate; Uncorrected includes no covariates. There are 142 genes significant in the same window in both models, and one shared gene significant only in different windows.'),
    dict(Item='Window definition', Details='Windows are symmetric TSS flanks of ±20, ±50, ±100, ±200, ±500 or ±1000 kb; full widths are twice the stated values. Separate columns list every significant window in each model. None means no FDR-significant window, not absence of a test.'),
    dict(Item='Regulatory correlations', Details='log10(species-median ATAC peak count) versus log10(MLS in years). The corrected model adds log10(body mass). These are counts of peaks overlapping the window, not integrated CPM restricted to peaks. Overview directions summarize significant windows only. Full effects/P/FDR/N for both models and every window are in Regulatory windows.'),
    dict(Item='Other assays retain correction', Details='Gene activity, liver transcriptome and liver proteome retain the body-mass-adjusted models from the original 216-gene table. These columns were not switched to the uncorrected model. This permits comparison of the regulatory selection groups against the same assay results.'),
    dict(Item='Direction symbols', Details='** means BH FDR < 0.05; * means nominal P < 0.05 only; unmarked Positive/Negative is not nominally significant. Not tested means no eligible fitted result. Not FDR-significant in a regulatory column means no window passes FDR in that model.'),
    dict(Item='Gene activity measures', Details='Unweighted integrated bigWig CPM·bp across TSS ±2 kb (4 kb total), or gene body plus 2 kb upstream promoter. Use log10 of species-median activity. Both adjust for log10(body mass); body also adjusts for log10(gene-body length). No pseudocount.'),
    dict(Item='Transcriptome', Details='Liver species median of source log.tpm, with no additional log transform, versus log10(MLS), adjusted for log10(body mass). BMR_Turk_Lu excluded; mouse liver supplement retained. ASM source common-name mismatch is retained from the existing species mapping.'),
    dict(Item='Proteome', Details='Liver species median of source quant.norm normalized on log.quant; no additional log transform. Correlation with log10(MLS) adjusts for log10(body mass).'),
    dict(Item='Species coverage', Details='Human, macaque and bowhead whale excluded; missing/lower-bound MLS excluded. Minimum five species per fit. Species vary by gene, window and assay and are not forced to match across assays. Exact IDs and N are retained in the detailed sheets.'),
    dict(Item='ATAC eligibility', Details='Choose the highest species-median locus separately by method before filtering, without fallback. Require positive values for log scales and >=50% of the intended span available within contig bounds. Internal gaps and read coverage are not measured by this rule.'),
    dict(Item='FDR families', Details='BH values are retained from all tested genes within each regulatory model/window and each genome-wide assay: TSS ±2 kb 20,244; gene body 20,200; RNA 12,804; protein 3,471. FDR is not recalculated over these 330 selected genes. Joint-across-window FDR is not the primary selection threshold.'),
    dict(Item='Aging evidence', Details=f'{supported} genes match at least one compiled evidence source. Yes means found in that collection; No is not proof that no aging study exists. GenAge includes human candidates and model homologs; CellAge includes manipulation and expression evidence; GWAS gene assignments are not necessarily causal.'),
    dict(Item='Experimental evidence', Details='Experiment flags the existing curated animal lifespan-extension collection, not all experiments in GenAge/CellAge. Human-homolog labels retain the tested organism and perturbation in Aging evidence. Activated fly Thor maps to EIF4EBP1; the finding does not imply that ordinary wild-type overexpression or higher human expression extends lifespan.'),
    dict(Item='Role definitions', Details='Roles come from source evidence, never project MLS correlation signs. Pro-longevity: loss/reduction shortens lifespan and/or overexpression extends it. Anti-longevity: loss/reduction extends lifespan and/or overexpression shortens it. Necessary-for-fitness and unannotated records are not forced into either role. https://genomics.senescence.info/help.html'),
    dict(Item='Cellular aging roles', Details='CellAge manipulation Induces means pro-senescence; Inhibits means anti-senescence. Expression signatures alone and GWAS associations do not establish causal aging direction. Cellular senescence is distinguished from whole-animal aging/lifespan. https://genomics.senescence.info/cells/index.html'),
    dict(Item='Mixed and unknown roles', Details='Conflicting directional evidence retains Mixed / context-dependent. Lifespan and cellular-senescence endpoints remain separate. Unknown means no directional evidence in the compiled collection.'),
    dict(Item='Symbols', Details='Original regulatory symbols are retained alongside HGNC-approved human symbols. Assays match exact source symbol first, then approved symbol, then unique HGNC alias/previous symbol. Ambiguous matches are not forced. Matching symbols does not establish one-to-one orthology.'),
    dict(Item='Interpretation', Details='Observational associations without phylogenetic, sequencing-depth, global peak-count or batch correction. Body-mass adjustment changes the association being estimated; neither selection group alone establishes causality. Source-organism or cell evidence does not prove human lifespan effects.'),
    dict(Item='Verification', Details=f'All original 216 gene rows retain the same evidence, roles and assay directions. {len(audit)} available assay fits were independently checked by full-model OLS. Genome-wide BH families, union membership, regulatory windows and source-file hashes were checked.'),
    dict(Item='Analysis sources', Details='Regulatory six-window analyses dated 22 October 2026 (corrected) and 23 October 2026 (uncorrected). Activity/omics genome-wide analysis dated 20 October 2026, mass_adjusted/all_available. Existing 216-gene integration and compiled aging/longevity reference dated 22 October 2026; database/study URLs are retained in Aging evidence.'),
]
data = dict(overview=overview, details=details, regulatory=regulatory, evidence=evidence, roles=roles, notes=notes,
            categories=categories, supported_genes=supported, source_counts=dict(Counter(k for r in overview for k in ['Experiment','GenAge','CellAge','GWAS'] if r[k]=='Yes')),
            role_summary=dict(directional_genes=sum(r['Directional_evidence'] for r in roles), longevity=dict(Counter(r['Longevity_role'] for r in roles)), cellular_senescence=dict(Counter(r['Cellular_senescence_role'] for r in roles))))
for name, rows in [('330_genes_overview.tsv',overview),('correlation_details.tsv',details),('regulatory_windows.tsv',regulatory),('aging_evidence.tsv',evidence),('aging_role_details.tsv',roles),('methods.tsv',notes)]:
    write(name,rows)
(OUT/'results.json').write_text(json.dumps(data,indent=2))
input_paths += [REG/'significant_genes_all_windows.tsv',RAW/'significant_genes_all_windows.tsv',REF/'aging_longevity_genes.tsv',REF/'evidence_details.tsv',REF/'sources/hgnc_complete_set.tsv',PREV/'results.json']
(OUT/'validation.json').write_text(json.dumps(dict(genes=len(genes),categories=categories,original_216_preserved=True,
    OLS=audit, source_sha256={str(p):sha(p) for p in input_paths}),indent=2))
print(json.dumps(dict(categories=categories, supported_genes=supported, source_counts=data['source_counts'],role_summary=data['role_summary'],evidence_rows=len(evidence),regulatory_rows=len(regulatory),assay_rows=len(details),recovered_symbols=[r for r in details if r['Symbol_matching']=='Unique HGNC previous/alias symbol']),indent=2))
