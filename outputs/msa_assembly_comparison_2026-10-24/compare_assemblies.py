"""Audit current browser reference builds against published MSA build metadata."""
import csv
import gzip
import html
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SRC = OUT / 'sources'
MULTIZ_SOURCE = 'https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons470way'
CACTUS_SOURCE = 'https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/trackDb/human/hg38/cactus447way.html'
TOGA_SOURCE = 'https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv'
ZOONOMIA_SOURCE = 'https://karlssonlab.org/zoonomia/'

species = list(csv.DictReader((ROOT / 'outputs/cactus447_species_overlap_2026-10-24/species_comparison.tsv').open(), delimiter='\t'))
inventory = {r['id']: r for r in json.loads((ROOT / 'outputs/reference_inventory_20261016/inventory.json').read_text())['species']}
cactus_table = json.loads((SRC / 'cactus447_assemblies.json').read_text())
multiz_table = json.loads((SRC / 'multiz470_assemblies.json').read_text())
ncbi = {r['accession']: r for r in json.loads((SRC / 'msa_ncbi_metadata.json').read_text())['reports']}
toga = {r['Assembly name']: r for r in csv.DictReader((SRC / 'toga_assembly_sources.tsv').open(), delimiter='\t')}
local_ncbi = {}
for path in (ROOT / 'outputs/reference_inventory_20261016/ncbi').glob('*.json'):
    for record in json.loads(path.read_text()).get('reports', []):
        if record.get('accession'):
            local_ncbi[record['accession']] = record

same_cactus = {'Human', 'Mouse', 'BMR', 'Chinchilla_lanigera', 'Octodon_degus', 'Bramus_lutescens', 'Scalopus_aquaticus', 'Ondatra_zibethicus'}
same_multiz = same_cactus | {'DMR', 'Molossus_molossus'}
partial_multiz = {'Sciurus_carolinensis'}
cactus_alias_accessions = {'Human': 'GCA_000001405.15', 'Mouse': 'GCA_000001635.2', 'DMR': 'GCF_000743615.1', 'Macaque': 'GCF_003339765.1', 'Chinchilla_lanigera': 'GCF_000276665.1'}

def record_for_acc(acc, code=''):
    record = ncbi[acc]
    return {'sequence_id': code, 'assembly': record['assembly_info']['assembly_name'], 'accession': acc,
            'source_url': f'https://www.ncbi.nlm.nih.gov/datasets/genome/{acc}/'}

def status(present, sid, same, partial=()):
    if not present:
        return 'Absent'
    if sid in same:
        return 'Matching build'
    if sid in partial:
        return 'Chromosomes match; scaffold version differs'
    return 'Different build / coordinates'

rows = []
evidence = []
for sp in species:
    sid = sp['browser_id']
    loc = inventory[sid]
    scientific = 'Ellobius lutescens' if sid == 'Bramus_lutescens' else sp['browser_scientific_name']
    # Inventory predates the activation of three species; current membership is
    # taken from the current 34-species catalog, not the inventory's status column.
    fasta = ROOT / loc['fasta']
    fai = Path(str(fasta) + '.fai')
    seqs = [line.split('\t') for line in fai.read_text().splitlines()]
    assert len(seqs) == loc['sequence_count']
    assert sum(int(s[1]) for s in seqs) == loc['total_bp']
    with (gzip.open(fasta, 'rt') if fasta.name.endswith('.gz') else fasta.open()) as handle:
        first_header = handle.readline().strip()
    evidence.append({'species': sid, 'fasta': str(fasta), 'first_header': first_header,
                     'sequence_count': len(seqs), 'total_bp': loc['total_bp'], 'inventory_dimensions_match': True})

    ct = [r for r in cactus_table if r['scientific'].split('/')[0] == scientific]
    assert bool(ct) == (sp['in_cactus447'] == 'Yes')
    assert len(ct) <= 1
    ca = []
    if ct:
        acc = re.findall(r'GC[AF]_\d+\.\d+', ct[0]['assembly'])
        if not acc:
            acc = [cactus_alias_accessions[sid]]
        ca = [record_for_acc(a, sp['cactus_sequence_ids']) for a in acc]
    mt = [r for r in multiz_table if r['cells'][3] == scientific or r['cells'][3].startswith(scientific + ' ')]
    assert bool(mt) == (sp['in_multiz470'] == 'Yes')
    ma = []
    if mt:
        for code in sp['multiz_assembly_ids'].split('; '):
            acc_source = 'GCA_000001405.15' if code == 'hg38' else toga[code]['NCBI accession / source']
            acc = re.findall(r'GC[AF]_\d+\.\d+', acc_source)
            if acc:
                ma.append(record_for_acc(acc[0], code))
            else:
                url = acc_source if acc_source.startswith('http') else 'https://www.dnazoo.org/assemblies/Castor_canadensis'
                label = 'Original bowhead scaffold assembly' if sid == 'Balaena_mysticetus' else 'DNA Zoo Hi-C assembly'
                ma.append({'sequence_id': code, 'assembly': label, 'accession': 'Not assigned in MSA source', 'source_url': url})

    cs = status(bool(ca), sid, same_cactus)
    ms = status(bool(ma), sid, same_multiz, partial_multiz)
    notes = []
    if sid in same_cactus or sid in same_multiz:
        notes.append('Matching nuclear assembly build from source metadata; full FASTA checksums were not compared. Confirm chromosome/scaffold aliases and sequence availability before reusing BED coordinates.')
    if sid == 'Human':
        notes.append('Local FASTA is the hg38 primary set. GENCODE v50/p14 describes the annotation release and does not make this a T2T human reference.')
    if sid == 'Mouse':
        notes.append('Use the mm10 row in Multiz470; its additional mm39 row is a different assembly. Cactus uses mm10/GRCm38.')
    if sid == 'DMR':
        notes.append('Multiz uses DMR_v1.0_HiC, matching the local nuclear build. Cactus uses older DMR_v1.0/fukDam1, without the Hi-C scaffolding.')
    if sid == 'Bramus_lutescens':
        notes.append('Bramus lutescens and Ellobius lutescens refer to the same species; both MSAs use ASM168507v1.')
    if sid == 'Balaena_mysticetus':
        notes.append('Local manifest explicitly names HLbalMyst1A and HAP1_SUPER sequences. This is different from Multiz HLbalMys1. Local public accession unresolved; Cactus lacks bowhead whale.')
    if sid == 'Myocastor_coypus':
        notes.append('Local Hi-C reference: 2,632,872,920 bp / 2,956 sequences; MSA MyoCoy_v1_BIUU: 2,903,175,784 bp. Exact public accession of local reference unresolved; do not treat these coordinate systems as identical.')
    if sid == 'Cricetulus_griseus':
        notes.append('Local GCF_003668045.3 pairs with GCA_003668045.2 (PICRH). Multiz uses GCA_003668045.1 (PICR); these are different scaffolding versions. Cactus uses another assembly, CHOK1S_HZDv1.')
    if sid == 'Sciurus_carolinensis':
        notes.append('Local GCF_902686445.1 pairs with GCA_902686445.2/mSciCar1.2; Multiz uses GCA_902686445.1/mSciCar1.1. All 21 assembled chromosome GenBank accessions and lengths match. All 752 sequence names/lengths match, but the other 731 scaffolds have different WGS accession versions and were not base-compared.')
    if sid == 'Molossus_molossus':
        notes.append('Hiller TOGA metadata resolves HLmolMol2 (Bat1K published) to GCA_014108415.1, paired to local GCF_014108415.1.')
    if sid == 'Tamias_striatus':
        notes.append('Local ptg contig assembly has no verified public accession; species absent from both alignments.')
    if sid == 'NMR':
        notes.append('Local chr-named mHetGlaV3 is not HetGla_female_1.0/hetGla2.')
    if sid == 'Peromyscus_maniculatus':
        notes.append('Local organism is P. maniculatus bairdii; neither old Pman_1.0 nor HU_Pman_2.1 is the local HU_Pman_BW_mat_3.1 assembly.')

    row = {'browser_id': sid, 'common_name': sp['common_name'], 'scientific_name': sp['browser_scientific_name'],
           'local_assembly': loc['assembly'], 'local_accession': loc['accession'] or 'Unresolved',
           'local_paired_genbank': local_ncbi.get(loc['accession'], {}).get('assembly_info', {}).get('paired_assembly', {}).get('accession', ''),
           'local_fasta': str(fasta), 'local_assembly_source': loc['assembly_url'],
           'cactus447_assembly': '; '.join(r['assembly'] for r in ca) or 'Absent',
           'cactus447_accession': '; '.join(r['accession'] for r in ca), 'cactus447_status': cs,
           'cactus447_sequence_id': sp['cactus_sequence_ids'],
           'multiz470_assembly': '; '.join(f"{r['sequence_id']}: {r['assembly']}" for r in ma) or 'Absent',
           'multiz470_accession': '; '.join(r['accession'] for r in ma), 'multiz470_status': ms,
           'notes': ' '.join(notes),
           'cactus447_assembly_links': '; '.join(r['source_url'] for r in ca),
           'multiz470_assembly_links': '; '.join(r['source_url'] for r in ma),
           'cactus447_table_source': CACTUS_SOURCE, 'multiz470_table_source': MULTIZ_SOURCE,
           'hiller_assembly_alias_source': TOGA_SOURCE,
           'cactus_dmr_accession_source': ZOONOMIA_SOURCE if sid == 'DMR' else ''}
    rows.append(row)

assert len(rows) == 34
assert len(same_cactus) == 8 and len(same_multiz) == 10
summary = {key: dict(Counter(r[key + '_status'] for r in rows)) for key in ['cactus447', 'multiz470']}
summary['scope'] = 'Current 34 browser species; assembly-build/coordinate comparison, not a full sequence identity or locus coverage assessment.'
summary['local_reference_validation'] = 'All 34 local FASTA first headers were inspected; current FASTA index counts and total lengths agree with the existing reference inventory.'
with (OUT / 'assembly_comparison_all_34_species.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
(OUT / 'assembly_comparison.json').write_text(json.dumps(rows, indent=2) + '\n')
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
(OUT / 'local_fasta_verification.json').write_text(json.dumps(evidence, indent=2) + '\n')

def esc(value):
    return html.escape(str(value), quote=True)

def badge(value):
    css = 'same' if value == 'Matching build' else 'partial' if value.startswith('Chromosomes') else 'absent' if value == 'Absent' else 'different'
    return f'<span class="badge {css}">{esc(value)}</span>'

def assembly_cell(name, acc, links):
    accs = acc.split('; ') if acc else []
    urls = links.split('; ') if links else []
    refs = ' · '.join(f'<a href="{esc(u)}">{esc(a)}</a>' for a, u in zip(accs, urls))
    return f'<strong>{esc(name)}</strong><div class="sub">{refs}</div>'

body = []
for r in rows:
    body.append('<tr>' +
                f'<td><strong>{esc(r["common_name"])}</strong><div class="sub"><i>{esc(r["scientific_name"])}</i></div></td>' +
                '<td>' + assembly_cell(r['local_assembly'], r['local_accession'], r['local_assembly_source']) + '</td>' +
                '<td>' + assembly_cell(r['cactus447_assembly'], r['cactus447_accession'], r['cactus447_assembly_links']) + badge(r['cactus447_status']) + '</td>' +
                '<td>' + assembly_cell(r['multiz470_assembly'], r['multiz470_accession'], r['multiz470_assembly_links']) + badge(r['multiz470_status']) + '</td>' +
                f'<td>{esc(r["notes"])}</td></tr>')
document = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Browser versus MSA genome assemblies</title>
<style>body{font:14px/1.45 system-ui,sans-serif;margin:32px;color:#183047;background:#fff}h1{font-size:26px;margin-bottom:8px}.intro{max-width:1100px}.summary{display:flex;gap:24px;margin:24px 0}.summary>div{background:#f1f6fa;padding:18px;border:1px solid #cbd9e4;border-radius:6px;flex:1}.summary b{font-size:20px}table{border-collapse:collapse;width:100%;margin-top:16px;table-layout:fixed}th,td{border:1px solid #d7e1e8;padding:12px;vertical-align:top;overflow-wrap:anywhere}th{background:#193e58;color:#fff;text-align:left;position:sticky;top:0}tbody tr:nth-child(even){background:#f7fafc}a{color:#075985}.sub{font-size:12px;margin-top:6px;color:#526476}.badge{display:inline-block;margin-top:8px;padding:3px 7px;border-radius:3px;font-size:12px}.same{background:#dff2e7;color:#125b36}.different{background:#fff0dc;color:#754400}.partial{background:#e9e5f8;color:#493477}.absent{background:#eaedf0;color:#586472}input{padding:10px;width:340px;border:1px solid #849aac;border-radius:4px}footer{margin:24px 0;color:#526476}</style></head><body>
<h1>Browser versus MSA genome assemblies</h1><p class="intro">34 current browser species · checked 24 October 2026. A species being present in an MSA does not mean its assembly matches your ATAC reference. “Matching build” means the same nuclear reference build based on assembly metadata, not a full FASTA checksum comparison. Chromosome aliases still need checking.</p>
<div class="summary"><div><b>Cactus447</b><br>23 species present · 8 matching builds · 15 different builds/coordinate systems · 11 absent</div><div><b>Multiz470</b><br>29 species present · 10 matching builds · 1 chromosome-level partial match · 18 different builds/coordinate systems · 5 absent</div></div>
<p>Eastern gray squirrel: all 21 assembled chromosome accessions match; other scaffold accession versions differ. Nutria and bowhead: the exact public accession of the local reference remains unresolved.</p>
<label for="filter">Find species, assembly or status: </label><input id="filter" placeholder="e.g. mouse, Different, GRCr8">
<table><colgroup><col style="width:14%"><col style="width:19%"><col style="width:19%"><col style="width:21%"><col style="width:27%"></colgroup><thead><tr><th>Species</th><th>Your browser reference</th><th>Cactus447 reference</th><th>Multiz470 reference</th><th>Interpretation</th></tr></thead><tbody>''' + ''.join(body) + '''</tbody></table>
<footer>Sources: <a href="''' + CACTUS_SOURCE + '''">UCSC Cactus447 assembly table</a> · <a href="''' + MULTIZ_SOURCE + '''">UCSC Multiz470 assembly table</a> · <a href="''' + TOGA_SOURCE + '''">Hiller assembly metadata</a> · NCBI assembly/sequence reports. Local references were checked against FASTA headers, index dimensions, the browser species catalog, and the existing reference inventory. No browser settings or data were changed.</footer>
<script>document.getElementById('filter').addEventListener('input',e=>{const q=e.target.value.toLowerCase();document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q))})</script></body></html>'''
(OUT / 'assembly_comparison.html').write_text(document)
print(json.dumps(summary, indent=2))
for r in rows:
    print(r['common_name'], '|', r['local_assembly'], '|', r['cactus447_assembly'], '|', r['cactus447_status'], '|', r['multiz470_assembly'], '|', r['multiz470_status'])
