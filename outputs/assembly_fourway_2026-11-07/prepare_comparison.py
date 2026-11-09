"""Combine the audited local/MSA inventory with the 2026-11-07 TOGA2 catalogue."""
import csv
import json
import re
import shutil
from collections import Counter
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
OLD = ROOT / 'outputs/msa_assembly_comparison_2026-10-24'
rows = list(csv.DictReader((OLD / 'assembly_comparison_all_34_species.tsv').open(), delimiter='\t'))
catalog = list(csv.DictReader((OUT / 'sources/toga2_assemblies.tsv').open(), delimiter='\t'))
inventory = {r['id']: r for r in json.loads((ROOT / 'outputs/reference_inventory_20261016/inventory.json').read_text())['species']}
by_code = {r['Assembly name']: r for r in catalog}
selected = {
    'Human': 'hg38', 'NMR': 'HLhetGla5', 'BMR': 'nanGal1', 'DMR': 'HLfukDam2',
    'Mouse': 'mm10', 'Rat': 'HLratNor8', 'Macaque': 'HLmacMula11', 'Rabbit': 'HLsylFlo1',
    'ASM': 'HLacoCah2', 'Hydrochoerus_hydrochaeris': 'HLhydHyd2',
    'Eonycteris_spelaea': 'HLeonSpe2A', 'Myocastor_coypus': 'HLmyoCoy2',
    'Balaena_mysticetus': 'HLbalMys1', 'Artibeus_jamaicensis': 'HLartJam3',
    'Castor_canadensis': 'HLcasCana4A', 'Cricetulus_griseus': 'HLcriGri3',
    'Chinchilla_lanigera': 'chiLan1', 'Octodon_degus': 'octDeg1',
    'Bramus_lutescens': 'HLellLut1', 'Eptesicus_fuscus': 'HLeptFus3',
    'Sciurus_niger': 'HLsciNig1', 'Meriones_unguiculatus': 'HLmerUng2',
    'Mesocricetus_auratus': 'HLmesAur3', 'Molossus_molossus': 'HLmolMol2',
    'Erethizon_dorsatum': 'HLereDor3', 'Tamiasciurus_hudsonicus': 'HLtamHuds1',
    'Marmota_monax': 'HLmarMon3', 'Peromyscus_maniculatus': 'HLperManiBai4',
    'Microtus_pennsylvanicus': 'HLmicPen1A', 'Cavia_porcellus': 'HLcavPor4',
    'Sciurus_carolinensis': 'HLsciCar1', 'Scalopus_aquaticus': 'HLscaAqu1',
    'Ondatra_zibethicus': 'HLondZib1', 'Tamias_striatus': None,
}
special = {
    'Human': 'Primary hg38 reference build matches. TOGA2 lists GRCh38.p12; the local FASTA is the primary hg38 set, so patch/alternate content is not asserted identical.',
    'Mouse': 'Use the mm10 entry in Multiz. Multiz also contains mm39, which is a different build. TOGA2 selection here is mm10.',
    'NMR': 'Local mHetGlaV3 differs from TOGA2 HLhetGla5 and the older hetGla2 build in both MSAs. TOGA2 also lists hetGla2 and DNA Zoo HLhetGla3.',
    'DMR': 'Multiz and TOGA2 match the local Hi-C build. Cactus uses the older DMR_v1.0/fukDam1 build.',
    'ASM': 'All 129 TOGA2 scaffold lengths match the local FASTA after mapping JAULSH010000NNN.1 to scaffold_NNN. Both MSAs use the older AcoCah_v1_BIUU.',
    'Hydrochoerus_hydrochaeris': 'All 343,878 sequence names and lengths match TOGA2 HLhydHyd2. Both MSAs use HydHyd_v1_BIUU, a different assembly.',
    'Myocastor_coypus': 'All 2,956 sequence names and lengths match TOGA2 HLmyoCoy2 (2,632,872,920 bp). Both MSAs use the older MyoCoy_v1_BIUU.',
    'Balaena_mysticetus': 'Local HLbalMyst1A uses HAP1_SUPER sequences and differs from the older HLbalMys1 in TOGA2 and Multiz. Absent from Cactus.',
    'Cricetulus_griseus': 'Local GCF_003668045.3 pairs with GenBank GCA_003668045.2. TOGA2 and Multiz use .1. Cactus uses a different CHO assembly.',
    'Sciurus_carolinensis': 'Local mSciCar1.2 differs in version from TOGA2/Multiz mSciCar1.1. All 21 assembled chromosome accessions and lengths match; other scaffold versions differ.',
    'Bramus_lutescens': 'Bramus lutescens and Ellobius lutescens are taxonomic aliases. All four datasets use ASM168507v1.',
    'Peromyscus_maniculatus': 'Local/TOGA2 matching accession is P. maniculatus bairdii. Do not substitute Peromyscus keeni or the older MSA assemblies.',
    'Tamias_striatus': 'Local ptg contig assembly has no verified public accession. No entry for this species in TOGA2, Multiz470 or Cactus447.',
}
size_keys = {'ASM': 'asm', 'Hydrochoerus_hydrochaeris': 'capy', 'Myocastor_coypus': 'nutria'}
for sid, key in size_keys.items():
    shutil.copy2(Path('/tmp') / f'nsc_toga2_{key}.sizes', OUT / 'sources' / f'{selected[sid]}.chrom.sizes')

output = []
for r in rows:
    sid = r['browser_id']
    loc = inventory[sid]
    fai = Path(str(ROOT / loc['fasta']) + '.fai')
    local_sizes = {p[0]: int(p[1]) for p in (line.split('\t') for line in fai.read_text().splitlines())}
    assert len(local_sizes) == loc['sequence_count']
    assert sum(local_sizes.values()) == loc['total_bp']
    c = by_code[selected[sid]] if selected[sid] else None
    support_url = ''
    if c is None:
        assert not any(x['Species'] == r['scientific_name'] for x in catalog)
        status, basis = 'Absent', 'Species absent from the public TOGA2 catalogue'
    elif sid in size_keys:
        remote = {p[0]: int(p[1]) for p in (line.split() for line in (OUT / 'sources' / f'{selected[sid]}.chrom.sizes').read_text().splitlines())}
        mapped = {f'scaffold_{int(k[8:].split(".")[0])}': v for k, v in local_sizes.items()} if sid == 'ASM' else local_sizes
        assert mapped == remote
        status = 'Matching build'
        basis = f'All {len(remote):,} sequence names/lengths match' + (' after documented scaffold alias conversion' if sid == 'ASM' else '')
        support_url = 'https://genome.senckenberg.de/download/TOGA2/Genomes/Mammalia/' + c['Directory Name'] + '/chrom.sizes'
    elif c['NCBI accession'] == r['local_accession']:
        status, basis = 'Matching build', 'Exact assembly accession and version'
    elif r['local_paired_genbank'] and c['NCBI accession'] == r['local_paired_genbank']:
        status, basis = 'Matching build', 'Recorded RefSeq/GenBank paired assemblies'
    elif sid in {'Human', 'Mouse'}:
        status, basis = 'Matching build', 'Same named nuclear reference build; complete FASTA contents not compared'
    elif sid == 'Sciurus_carolinensis':
        status, basis = 'Different version; chromosomes match', 'Different assembly version, but all 21 assembled chromosome accessions/lengths match'
    else:
        status, basis = 'Different build/version', 'No local build match in catalogue; selected related assembly differs'
    for source in ['multiz470', 'cactus447']:
        r[source+'_status'] = r[source+'_status'].replace('Different build / coordinates', 'Different build/version').replace('Chromosomes match; scaffold version differs', 'Different version; chromosomes match')
    acc = c['NCBI accession'] if c and re.fullmatch(r'GC[AF]_\d+\.\d+', c['NCBI accession']) else 'Not assigned in catalogue' if c else ''
    r.update({
        'toga2_assembly_id': selected[sid] or 'Absent',
        'toga2_assembly': (c['NCBI assembly name'] or {'HLacoCah2':'Acomys chromosome-scale assembly', 'HLhydHyd2':'DNA Zoo Hi-C assembly','HLmyoCoy2':'DNA Zoo Hi-C assembly','HLbalMys1':'Original bowhead scaffold assembly'}.get(selected[sid],selected[sid])) if c else 'Absent',
        'toga2_accession': acc,
        'toga2_status': status,
        'toga2_match_basis': basis,
        'toga2_genome_url': 'https://genome.senckenberg.de/download/TOGA2/Genomes/Mammalia/' + c['Directory Name'] + '/' if c else '',
        'toga2_catalogue_url': 'https://genome.senckenberg.de/download/TOGA2/assemblies_and_species.tsv',
        'support_url': support_url,
        'coordinate_note': special.get(sid, 'Compare chromosome/scaffold names before transferring coordinates. Matching build metadata does not establish full sequence identity.'),
        'checked_date': '2026-11-07',
    })
    if sid == 'Bramus_lutescens':
        r['common_name'] = 'Transcaucasian mole vole'
    output.append(r)

assert len(output) == len({r['browser_id'] for r in output}) == 34
summary = {}
for source in ['toga2','multiz470','cactus447']:
    counts = Counter(r[source+'_status'] for r in output)
    summary[source] = {'present':34-counts['Absent'], 'matching':counts['Matching build'],
        'different':counts['Different build/version'] + counts['Different version; chromosomes match'],
        'chromosome_match_subset':counts['Different version; chromosomes match'], 'absent':counts['Absent']}
assert summary['toga2'] == dict(present=33,matching=28,different=5,chromosome_match_subset=1,absent=1)
assert summary['multiz470'] == dict(present=29,matching=10,different=19,chromosome_match_subset=1,absent=5)
assert summary['cactus447'] == dict(present=23,matching=8,different=15,chromosome_match_subset=0,absent=11)
for x in summary.values(): assert x['matching']+x['different']+x['absent']==34
all_match = [r['common_name'] for r in output if all(r[k+'_status']=='Matching build' for k in summary)]
assert len(all_match)==8
(OUT/'comparison.json').write_text(json.dumps({'rows':output,'summary':summary,'matching_all_four':all_match},indent=2)+'\n')
with (OUT/'Assembly_comparison_34_species.tsv').open('w',newline='') as f:
    fields=['common_name','scientific_name','local_assembly','local_accession','local_paired_genbank',
            'toga2_assembly_id','toga2_assembly','toga2_accession','toga2_status','toga2_match_basis',
            'multiz470_assembly','multiz470_accession','multiz470_status',
            'cactus447_assembly','cactus447_accession','cactus447_status','coordinate_note',
            'local_assembly_source','toga2_genome_url','toga2_catalogue_url',
            'multiz470_table_source','multiz470_assembly_links','cactus447_table_source','cactus447_assembly_links','support_url','checked_date']
    writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); writer.writeheader();writer.writerows(output)
print(json.dumps({'summary':summary,'matching_all_four':all_match},indent=2))
