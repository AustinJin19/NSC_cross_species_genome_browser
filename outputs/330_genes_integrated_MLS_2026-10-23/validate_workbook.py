"""Read-only verification of the saved Excel output against compiled results."""
import json
import math
from collections import Counter
from pathlib import Path
import openpyxl

OUT = Path(__file__).resolve().parent
d = json.loads((OUT/'results.json').read_text())
p = OUT/'330_genes_corrected_and_uncorrected_MLS.xlsx'
wb = openpyxl.load_workbook(p, data_only=True)
ws = wb['330 genes']
keys = ['Gene','Human_symbol','Regulatory_significance','Corrected_windows_kb','Uncorrected_windows_kb',
        'Experiment','GenAge','CellAge','GWAS','Reported_aging_longevity_role',
        'Regulatory_corrected','Regulatory_uncorrected','tss2kb','body','rna','protein']
assert ws.max_row==336 and ws.max_column==16
for index, record in enumerate(d['overview'],7):
    assert [ws.cell(index,j+1).value for j in range(len(keys))] == [record[k] for k in keys], record['Gene']
assert Counter(ws.cell(i,3).value for i in range(7,337))==Counter(d['categories'])
assert ws.freeze_panes=='C7'
assert '143' in ws['A3'].value and '73' in ws['A3'].value and '114' in ws['A3'].value
assert ws.tables['CombinedGeneOverview'].ref=='A6:P336'
assert ws.tables['CombinedGeneOverview'].autoFilter.ref=='A6:P336'
for sheet,key in [('Regulatory windows','regulatory'),('Correlation details','details'),('Aging evidence','evidence'),('Role details','roles'),('Methods','notes')]:
    assert wb[sheet].max_row==len(d[key])+1, sheet

def equal(a,b):
    if b is None or b=='':
        return a is None or a==''
    if isinstance(b,(int,float)) and not isinstance(b,bool):
        return isinstance(a,(int,float)) and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)
    return a==b

reg_keys=['Gene','Human_symbol','Model','Body_mass_correction','Window_flank_kb','Full_width_kb','Direction',
          'Correlation_r','P','BH_FDR_window','Significant_FDR05','N_species','DF','MLS_coefficient','FDR_family','Species_IDs']
assay_keys=['Gene','Human_symbol','Method','Body_mass_correction','Covariates','Direction','Partial_r','P','BH_FDR',
            'N_species','DF','MLS_coefficient','Status','Source_gene_symbol','Symbol_matching','FDR_family','Species_IDs']
for name,records,cols in [('Regulatory windows',d['regulatory'],reg_keys),('Correlation details',d['details'],assay_keys)]:
    for i,r in enumerate(records,2):
        for j,k in enumerate(cols,1):
            assert equal(wb[name].cell(i,j).value,r[k]),(name,i,k)
for i,r in enumerate(d['evidence'],2):
    for j,k in enumerate(r,1):
        expected = float(r[k]) if k=='P value' and r[k] else r[k]
        assert equal(wb['Aging evidence'].cell(i,j).value,expected),(i,k)
role_keys=['Gene','Human_symbol','Reported_role','Longevity_role','Cellular_senescence_role','Directional_evidence','Evidence_basis','Interpretation','Source_URLs']
for i,r in enumerate(d['roles'],2):
    for j,k in enumerate(role_keys,1):
        expected=('Yes' if r[k] else 'No') if k=='Directional_evidence' else r[k]
        assert equal(wb['Role details'].cell(i,j).value,expected),(i,k)
errors=[(s.title,c.coordinate,c.value) for s in wb for row in s for c in row if c.data_type=='e']
assert not errors,errors
significant=Counter((r['Model'],r['Window_flank_kb']) for r in d['regulatory'] if r['Significant_FDR05']=='Yes')
assert [significant['Corrected',w] for w in [20,50,100,200,500,1000]]==[4,6,20,78,71,121]
assert [significant['Uncorrected',w] for w in [20,50,100,200,500,1000]]==[6,5,22,29,138,168]
report=dict(gene_rows=330,overview_cells_verified=330*16,regulatory_rows_verified=3960,assay_rows_verified=1320,
            evidence_rows_verified=len(d['evidence']),role_rows_verified=330,categories=d['categories'],
            filters=True,freeze_panes=ws.freeze_panes,summary_formula_cached=ws['A3'].value,formula_errors=errors)
(OUT/'workbook_validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
