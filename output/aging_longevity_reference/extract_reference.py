from pathlib import Path
import csv, json
import openpyxl
ROOT=Path(__file__).resolve().parent
wb=openpyxl.load_workbook(ROOT/'sources/Tyshkovskiy2026_Table2.xlsx',read_only=True,data_only=True)
url='https://www.nature.com/articles/s41586-026-10542-3'
def write(name,rows,fields):
    with (ROOT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
allrows=[]; bysheet={}; counts=[]
for ws in wb:
    it=ws.iter_rows(values_only=True);headers=next(it);rows=[]
    for vals in it:
        d=dict(zip(headers,vals));symbol=d.get('Gene.symbol');p=d.get('P.Adjusted');s=d.get('Slope')
        if not symbol or not isinstance(p,(int,float)) or not isinstance(s,(int,float)):continue
        r={'gene_symbol_source':symbol,'symbol_uppercase_lookup':str(symbol).upper(),'entrez_id_source':d['Entrez.ID'],'signature':ws.title,'slope':s,'direction':'positive' if s>0 else 'negative' if s<0 else 'zero','p_value':d['P.Value'],'adjusted_p_value':p,'source':url}
        rows.append(r)
    bysheet[ws.title]=rows
    sig=[r for r in rows if r['adjusted_p_value']<.05];allrows+=sig
    counts.append({'signature':ws.title,'tested_genes':len(rows),'significant_genes':len(sig),'positive':sum(r['slope']>0 for r in sig),'negative':sum(r['slope']<0 for r in sig)})
fields=list(allrows[0]);write('2026_all_significant_associations.tsv',allrows,fields)
write('2026_signature_counts.tsv',counts,list(counts[0]))
for prefix,name in [('(A)','2026_mouse_liver_aging'),('(G)','2026_mouse_liver_lifespan_age_adjusted'),('(N)','2026_rodent_lifespan_age_adjusted'),('(E)','2026_mouse_liver_mortality_age_adjusted')]:
    rows=[r for k,v in bysheet.items() if k.startswith(prefix) for r in v if r['adjusted_p_value']<.05]
    write(name+'.tsv',rows,fields)
    for direction in ['positive','negative']:
        (ROOT/(name+'_'+direction+'.txt')).write_text('\n'.join(sorted({r['symbol_uppercase_lookup'] for r in rows if r['direction']==direction}))+'\n')
maps=[{r['gene_symbol_source']:r for r in v} for k,v in bysheet.items() if k.startswith(('(O)','(P)','(Q)','(R)'))]
conserved=[]
for gene in sorted(set.intersection(*(set(m) for m in maps))):
    rr=[m[gene] for m in maps]
    if all(r['adjusted_p_value']<.05 for r in rr) and len({r['direction'] for r in rr})==1:
        conserved.append({'gene_symbol_source':gene,'symbol_uppercase_lookup':gene.upper(),'direction_with_age':rr[0]['direction'],**{label+'_adjusted_p':r['adjusted_p_value'] for label,r in zip(['mouse','rat','macaque','human'],rr)}})
write('2026_conserved_aging_strict_four_species.tsv',conserved,list(conserved[0]) if conserved else ['gene_symbol_source'])
print(json.dumps({'signature_counts':counts,'strict_four_species_count':len(conserved),'strict_examples':conserved[:10]},indent=2))
