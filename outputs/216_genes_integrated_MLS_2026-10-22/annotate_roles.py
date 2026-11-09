"""Source-based longevity/senescence roles; never infer roles from MLS correlations."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
data = json.loads((OUT / 'results.json').read_text())
evidence_by_gene = defaultdict(list)
for e in data['evidence']:
    resource, typ, direction = e['Resource'], e['Evidence type'], e['Direction in source']
    role, scope = 'Unclassified', 'No directional perturbation evidence in this record'
    if resource == 'Experiment' and typ == 'Animal lifespan extension':
        if direction == 'Increase':
            role = 'Pro-longevity'
        elif direction == 'Decrease':
            role = 'Anti-longevity'
        scope = 'Lifespan extension intervention in ' + e['Tested organism']
    elif resource == 'GenAge' and direction in ('Pro-Longevity', 'Anti-Longevity'):
        role = 'Pro-longevity' if direction == 'Pro-Longevity' else 'Anti-longevity'
        scope = 'Source-organism label: ' + e['Tested organism'] + '; human homolog inference only'
    elif resource == 'CellAge' and typ == 'Senescence manipulation':
        if direction == 'Induces':
            role = 'Pro-senescence'
        elif direction == 'Inhibits':
            role = 'Anti-senescence'
        scope = 'Human cellular senescence; not a demonstrated whole-animal aging/lifespan role'
    elif resource == 'GWAS':
        scope = 'Variant association; gene functional direction not established'
    elif typ == 'Senescence expression':
        scope = 'Expression association; pro-/anti-aging causality not established'
    elif direction == 'Necessary for fitness':
        scope = 'Fitness evidence; not reclassified as pro-longevity'
    e['Interpreted role'] = role
    e['Role scope'] = scope
    evidence_by_gene[e['Human gene']].append(e)

role_rows = []
for row in data['overview']:
    evidence = evidence_by_gene[row['Human_symbol']]
    life = {e['Interpreted role'] for e in evidence if e['Interpreted role'] in ('Pro-longevity', 'Anti-longevity')}
    cell = {e['Interpreted role'] for e in evidence if e['Interpreted role'] in ('Pro-senescence', 'Anti-senescence')}
    life_role = 'Mixed / context-dependent' if len(life) > 1 else next(iter(life)) if life else 'Unknown'
    cell_role = 'Mixed / context-dependent' if len(cell) > 1 else next(iter(cell)) if cell else 'Unknown'
    pieces = []
    if len(life) > 1:
        pieces.append('Mixed longevity roles (model homologs)')
    elif life:
        label = next(iter(life))
        scope = 'mouse experiment + model evidence' if row['Experiment'] == 'Yes' else 'model homolog evidence'
        pieces.append(f'{label} ({scope})')
    if len(cell) > 1:
        pieces.append('Mixed senescence roles (cells)')
    elif cell:
        pieces.append(next(iter(cell)) + ' (cellular aging only)')
    if not pieces:
        pieces = ['Unknown / not established' if evidence else 'Unknown (not in compiled reference)']
    role = '\n'.join(pieces)
    row['Reported_aging_longevity_role'] = role
    basis = '; '.join(sorted({f"{e['Resource']}: {e['Original gene']} ({e['Tested organism']}), {e['Direction in source']}"
                              for e in evidence if e['Interpreted role'] != 'Unclassified'}))
    if not basis and evidence:
        basis = 'Existing association, human-candidate or fitness evidence does not establish a pro-/anti-aging role.'
    role_rows.append(dict(Gene=row['Gene'], Human_symbol=row['Human_symbol'],
        Reported_role=role, Longevity_role=life_role, Cellular_senescence_role=cell_role,
        Directional_evidence=bool(life or cell), Evidence_basis=basis,
        Interpretation='Roles describe the cited source context; homologs and cellular senescence are not established human lifespan effects.',
        Source_URLs=' ; '.join(sorted({e['Source URL'] for e in evidence if e['Interpreted role'] != 'Unclassified'}))))

notes = [
    {'Item':'Reported aging/longevity role', 'Details':'The added column is derived only from source evidence, never from the sign of the project MLS correlation. Pro-longevity, anti-longevity and mixed refer to the source organism/model. Cellular aging roles are labelled pro-senescence or anti-senescence, explicitly confined to cellular assays; these do not establish whole-animal pro-/anti-aging effects.'},
    {'Item':'GenAge role definitions', 'Details':'Pro-longevity: loss/reduction shortens lifespan and/or overexpression extends it. Anti-longevity: loss/reduction extends lifespan and/or overexpression shortens it. Necessary for fitness, Unannotated and human-candidate records are not forced into either class. Definitions: https://genomics.senescence.info/help.html'},
    {'Item':'CellAge role definitions', 'Details':'Manipulation records annotated Induces are called pro-senescence; Inhibits are anti-senescence. Senescence-expression up/down signatures alone are not assigned a causal role. Cellular senescence can have different consequences across tissues and disease contexts. Source scope: https://genomics.senescence.info/cells/index.html'},
    {'Item':'Unknown and mixed roles', 'Details':'GWAS/proximity associations, expression signatures and absent evidence remain Unknown unless another directional source exists. CERS4 and TRIP13 retain mixed pro-/anti-longevity model-homolog annotations. KDM4C retains its anti-longevity model label and anti-senescence cellular label separately; distinct endpoints are not collapsed into one sign.'},
]
items = {n['Item'] for n in notes}
data['notes'] = [n for n in data['notes'] if n['Item'] not in items] + notes
data['roles'] = role_rows
data['role_summary'] = dict(directional_genes=sum(r['Directional_evidence'] for r in role_rows),
    longevity=dict(Counter(r['Longevity_role'] for r in role_rows)),
    cellular_senescence=dict(Counter(r['Cellular_senescence_role'] for r in role_rows)))

def write(name, rows):
    with (OUT / name).open('w') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

write('216_genes_overview.tsv', data['overview'])
write('aging_evidence.tsv', data['evidence'])
write('aging_role_details.tsv', role_rows)
write('methods.tsv', data['notes'])
(OUT / 'results.json').write_text(json.dumps(data, indent=2))
print(json.dumps(data['role_summary'], indent=2))
