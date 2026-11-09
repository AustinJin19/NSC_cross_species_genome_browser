import csv, gzip, json, pathlib

root = pathlib.Path('/Users/zjin27/Downloads/CS_Genome_Browser')
out = root / 'outputs/reference_inventory_20261016'
rows = []
for folder in sorted((root/'reference_genomes').glob('**/*')):
    if not folder.is_dir():
        continue
    files = sorted(folder.iterdir())
    fas = [p for p in files if p.name.endswith(('.fa','.fa.gz'))]
    anns = [p for p in files if p.name.endswith(('.gtf.gz','.gtf')) or p.name == 'Acomys_cahirinus.gz']
    if not fas and not anns:
        continue
    r = {'folder':str(folder.relative_to(root)), 'fasta':[], 'annotations':[]}
    for p in fas:
        op = gzip.open if p.suffix == '.gz' else open
        with op(p,'rt') as f:
            header = f.readline().strip()
        fai = pathlib.Path(str(p)+'.fai')
        seqs = []
        if fai.exists():
            seqs = [line.strip().split('\t')[:2] for line in fai.open()]
        r['fasta'].append({'path':str(p.relative_to(root)), 'header':header, 'sequences':len(seqs), 'total_bp':sum(int(x[1]) for x in seqs), 'first_seqs':seqs[:4]})
    for p in anns:
        op = gzip.open if p.suffix == '.gz' else open
        comments, features = [], []
        with op(p,'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    comments.append(line)
                elif line:
                    features.append(line)
                    if len(features) >= 2:
                        break
        r['annotations'].append({'path':str(p.relative_to(root)), 'comments':comments, 'features':features})
    rows.append(r)
(out/'local_evidence.json').write_text(json.dumps(rows,indent=2))
for r in rows:
    print('\n'+r['folder'])
    for f in r['fasta']:
        print('FASTA:',f['header'],'|',f['sequences'],'seqs',f['total_bp'],'bp','|',f['first_seqs'][:2])
    for a in r['annotations']:
        print('GTF:',pathlib.Path(a['path']).name,'|','; '.join(a['comments']))
        print('FEATURE:',a['features'][0][:550] if a['features'] else '')
