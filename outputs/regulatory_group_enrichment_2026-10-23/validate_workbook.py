import json
import math
from collections import Counter
from pathlib import Path
import openpyxl

root=Path(__file__).resolve().parent
d=json.loads((root/'report_data.json').read_text())
file=root/'Regulatory_pathway_enrichment_comparison.xlsx'
wb=openpyxl.load_workbook(file,data_only=True)
formulas=openpyxl.load_workbook(file,data_only=False)
assert wb.sheetnames==['Comparison','Significant terms','Plot terms','Input genes','Methods']
assert [[wb['Comparison'].cell(r,c).value for c in range(2,6)] for r in range(7,10)]==[[143,143,18,2],[73,71,7,0],[114,114,0,0]]
for sheet,key in [('Significant terms','significant'),('Plot terms','comparison')]:
    ws=wb[sheet]
    assert ws.max_row==len(d[key])+1
    for i,r in enumerate(d[key],2):
        expected=[r['Group'],r['Library'],r['Term_ID'],r['Term'],r['BH_FDR'],r['Fold_enrichment'],r['Overlap'],
                  r['Query_size'],r['Gene_ratio'],r['Term_background'],r['Background_size'],r['Raw_P'],
                  r['Joint_BH_FDR_six_queries'],'Yes' if r['Significant_FDR05'] else 'No',r['Genes'] or None,r['Source_URL']]
        for c,value in enumerate(expected,1):
            found=ws.cell(i,c).value
            if isinstance(value,(int,float)):
                assert isinstance(found,(int,float)) and math.isclose(found,value,rel_tol=1e-12,abs_tol=1e-15),(sheet,i,c,found,value)
            else:assert found==value,(sheet,i,c,found,value)
genes=list(wb['Input genes'].iter_rows(min_row=2,values_only=True))
assert len(genes)==330 and len({x[1] for x in genes})==330
assert Counter(x[0] for x in genes)==Counter({'Shared':143,'Corrected only':73,'Uncorrected only':114})
assert {x[1] for x in genes if x[8]!='Mapped'}=={'MS4A20','ZFP865'}
for x,r in zip(genes,d['input_genes']):
    assert x[1]==r['Gene'] and x[2]==r['Human_symbol'] and x[3]==r['Submitted_ID']
for ws in wb:
    assert len(ws.tables)==1
    t=next(iter(ws.tables.values()))
    assert t.autoFilter is not None and t.autoFilter.ref==t.ref
    assert ws.freeze_panes is not None
    for row in ws:
        for cell in row:
            assert cell.data_type!='e',(ws.title,cell.coordinate,cell.value)
for row in range(7,10):
    for column in range(2,6):assert formulas['Comparison'].cell(row,column).data_type=='f'
out={'validated':True,'sheets':wb.sheetnames,'significant_term_rows':27,'comparison_rows':36,
     'input_genes':330,'effective_background':19196,'mapped_genes':328,
     'summary_formula_cells':12,'tables_with_filters':5,'formula_errors':0}
(root/'workbook_validation.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
