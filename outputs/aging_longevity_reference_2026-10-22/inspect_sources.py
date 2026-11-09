import csv, io, zipfile, re, json
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parent
src = root / 'sources'
z = zipfile.ZipFile(src / 'gwas_catalog_full.zip')
print('GWAS archive', [(x.filename, x.file_size) for x in z.infolist()])
name = next(n for n in z.namelist() if not n.endswith('/'))
reader = csv.DictReader(io.TextIOWrapper(z.open(name), encoding='utf-8-sig'), delimiter='\t')
print('GWAS columns', reader.fieldnames)
terms = Counter()
kept = []
for row in reader:
    trait = row['DISEASE/TRAIT']
    mapped = row.get('MAPPED_TRAIT', '')
    if re.search(r'\b(longevity|lifespan|life span|healthspan|health span|aging|ageing|epigenetic age|biological age|age at death)\b', trait + ' ' + mapped, re.I):
        kept.append(row)
        terms[(trait, mapped)] += 1
(src / 'gwas_aging_search.json').write_text(json.dumps(kept, indent=2))
(src / 'gwas_trait_candidates.tsv').write_text('Reported trait\tMapped trait\tAssociation rows\n' + ''.join(f'{t}\t{m}\t{n}\n' for (t,m),n in sorted(terms.items())))
print('Candidate associations', len(kept), 'trait combinations', len(terms))
for (t,m),n in sorted(terms.items()): print(n,t,'|',m)
