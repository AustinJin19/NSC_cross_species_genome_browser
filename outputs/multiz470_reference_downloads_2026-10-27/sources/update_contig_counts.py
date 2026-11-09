"""Add source-reported, version-specific contig counts to the reference catalog."""
import csv
import html
import json
from pathlib import Path
import re
import zipfile

OUT = Path(__file__).resolve().parents[1]
ROOT = OUT.parents[1]
NCBI_SNAPSHOT = ROOT / 'outputs/msa_assembly_comparison_2026-10-24/sources/msa_ncbi_metadata.json'
metadata = {r['accession']: r for r in json.loads(NCBI_SNAPSHOT.read_text())['reports']}
rows = json.loads((OUT / 'resolved_sources.json').read_text())
e = html.escape

for row in rows:
    accession = row['accession']
    if accession in metadata:
        record = metadata[accession]
        assert row['assembly'] == record['assembly_info']['assembly_name']
        row.update(
            contig_count=record['assembly_stats']['number_of_contigs'],
            contig_count_status='NCBI reported',
            contig_count_source='NCBI Datasets assembly_stats.number_of_contigs',
            contig_count_source_url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report?filters.assembly_version=all_assemblies',
            contig_count_notes='Assembly-wide count for the exact accession.version, from the 2026-10-24 NCBI snapshot. Not a count of FASTA records or scaffolds; the downloadable file was not recounted.',
        )
    elif row['msa_id'] == 'HLbalMys1':
        row.update(
            contig_count=113673,
            contig_count_status='Published source assembly',
            contig_count_source='Keane et al. (2015), Cell Reports, Table 1',
            contig_count_source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4536333/#tbl1',
            contig_count_notes='Published count for the original bowhead source assembly. Not independently recomputed from the archived HLbalMys1.2bit; contig/gap definitions can differ between sources.',
        )
    elif row['msa_id'] == 'HLcasCan3':
        row.update(
            contig_count=None,
            contig_count_status='Not verified',
            contig_count_source='Hiller/DNA Zoo archive; exact-assembly count not found in checked sources',
            contig_count_source_url='https://genome.senckenberg.de/download/TOGA/MammalianDNAZooAssemblies/',
            contig_count_notes='No verified contig count found for this exact custom Hi-C assembly. Neither scaffold count nor the count from a different beaver assembly is substituted.',
        )
    else:
        raise ValueError(f"Missing count provenance: {row['msa_id']}")

assert len(rows) == 33
assert sum(row['contig_count'] is not None for row in rows) == 32
(OUT / 'resolved_sources.json').write_text(json.dumps(rows, indent=2) + '\n')
stats_fields = ['contig_count', 'contig_count_status', 'contig_count_source', 'contig_count_source_url', 'contig_count_notes']
evidence = [{k: r[k] for k in ['msa_id', 'accession', 'assembly'] + stats_fields} for r in rows]
(OUT / 'sources/assembly_contig_stats.json').write_text(json.dumps(evidence, indent=2) + '\n')

manifest = OUT / 'multiz470_reference_downloads.tsv'
with manifest.open() as stream:
    columns = next(csv.reader(stream, delimiter='\t'))
columns = [column for column in columns if column not in stats_fields]
pos = columns.index('accession') + 1
columns[pos:pos] = stats_fields
with manifest.open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', extrasaction='ignore', lineterminator='\n')
    writer.writeheader()
    writer.writerows({k: ('NA' if k == 'contig_count' and v is None else v) for k, v in row.items()} for row in rows)

tr = []
for row in rows:
    if row['ftp_url'] != '-':
        ftp = f'<a href="{e(row["ftp_url"])}">{e(row["ftp_url"])}</a>'
    else:
        ftp = 'HTTPS archive only' if row['status'] == 'available' else 'FASTA unavailable'
    if row['status'] == 'available':
        direct = f'<a href="{e(row["download_url"])}">{e(row["filename"])}</a>'
    else:
        direct = f'<a href="{e(row["directory"])}">NCBI release directory</a>'
    count = f'{row["contig_count"]:,}' if row['contig_count'] is not None else 'NA'
    source_label = {'NCBI reported': 'NCBI', 'Published source assembly': 'Paper†', 'Not verified': 'Not verified'}[row['contig_count_status']]
    common = re.sub(r' \((?:mm10|hg38)\)$', '', row['common_name'])
    size = f'{row["bytes"]/1e9:.3f} GB' if row['bytes'] else 'unavailable'
    cells = [
        f'<td>{e(common)}<br><small><i>{e(row["scientific_name"])}</i></small></td>',
        f'<td><b>{e(row["msa_id"])}</b><br>{e(row["assembly"])}<br><small>{e(row["accession"])}</small></td>',
        f'<td class="contigs" title="{e(row["contig_count_notes"])}"><strong>{count}</strong><br><small><a href="{e(row["contig_count_source_url"])}">{source_label}</a></small></td>',
        f'<td class="url">{ftp}</td>',
        f'<td class="download">{direct}<br><small>{e(row["format"])} · {size}</small></td>',
        f'<td>{e(row["notes"])}</td>',
    ]
    tr.append('<tr>' + ''.join(cells) + '</tr>')
table = '<table><thead><tr><th>Species</th><th>Multiz ID / assembly</th><th>Contig count</th><th>FTP path</th><th>Verified HTTPS download</th><th>Notes</th></tr></thead><tbody>\n' + '\n'.join(tr) + '\n</tbody></table>'
note = '''<p id="contig-count-note"><b>Contig counts:</b> 31 entries use NCBI’s assembly-wide <code>number_of_contigs</code> for the exact accession version (metadata retrieved 24 October 2026). These describe contigs, not scaffold totals or FASTA record counts; source definitions and sequence sets may differ. <b>† Bowhead:</b> 113,673 is the count for the original source assembly in <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4536333/#tbl1">Keane et al. (2015), Table 1</a>; the archived 2bit file was not recounted. <b>Beaver:</b> NA means no verified count was found for the exact HLcasCan3 assembly. Each reported count links to its source.</p>'''
report = OUT / 'reference_downloads.html'
page = report.read_text()
page = re.sub(r'<p id="contig-count-note">.*?</p>', '', page, flags=re.S)
page = re.sub(r'<table>.*?</table>', note + '\n' + table, page, flags=re.S)
if '.contigs{' not in page:
    page = page.replace('</style>', '.contigs{min-width:100px;text-align:right;font-variant-numeric:tabular-nums}.download{overflow-wrap:anywhere;max-width:300px}</style>')
report.write_text(page)

validation_path = OUT / 'sources/validation.json'
validation = json.loads(validation_path.read_text())
validation['contig_counts'] = {'ncbi_reported': 31, 'published_source_assembly': 1, 'not_verified': ['HLcasCan3'], 'updated_on': '2026-10-27'}
validation_path.write_text(json.dumps(validation, indent=2) + '\n')
with zipfile.ZipFile(OUT / 'multiz470_reference_download_package.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
    for name in ['download_multiz470_references.sh', 'multiz470_reference_downloads.tsv', 'ftp_paths.txt', 'reference_downloads.html']:
        archive.write(OUT / name, arcname='multiz470_reference_downloads/' + name)
    archive.write(validation_path, arcname='multiz470_reference_downloads/validation.json')
    archive.write(OUT / 'sources/assembly_contig_stats.json', arcname='multiz470_reference_downloads/assembly_contig_stats.json')
print('Updated HTML, TSV, JSON and ZIP: 31 NCBI counts, 1 published source count, 1 explicitly unavailable.')
