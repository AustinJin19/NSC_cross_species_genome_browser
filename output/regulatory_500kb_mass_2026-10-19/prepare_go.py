import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent / 'outputs/regulatory-500kb-mass-2026-10-19'
GO = HERE / 'go'
GO.mkdir(exist_ok=True)
with (OUT/'all_tested.tsv').open() as f:
    background=[r['gene'] for r in csv.DictReader(f, delimiter='\t')]
assert len(background)==len(set(background))
manifest=[]
for threshold in ['FDR05','nominal_p05']:
    for direction in ['positive','negative']:
        name=f'{threshold}_{direction}'
        genes=(OUT/f'{name}_genes.txt').read_text().splitlines()
        assert len(genes)==len(set(genes)) and set(genes)<=set(background)
        manifest.append(dict(name=name, threshold=threshold, direction=direction,
                             input_genes=len(genes), submitted_background=len(background), run=bool(genes)))
        if genes:
            request=dict(organism='hsapiens', query=genes, background=background,
                         domain_scope='custom', sources=['GO:BP'], significance_threshold_method='fdr',
                         user_threshold=.05, all_results=True, no_evidences=False, ordered=False)
            (GO/f'{name}_request.json').write_text(json.dumps(request))
(GO/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest,indent=2))
