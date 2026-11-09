"""Select diverse finished fits for independent browser/R validation."""
import csv,gzip,json,sys
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
from analyze import METHODS

def main():
    cases=[];comparisons=[];rng=np.random.default_rng(27)
    for method in METHODS:
        file=OUT/'results'/method/'all_tested.tsv.gz'
        if not file.exists():continue
        with gzip.open(file,'rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
        with np.load(OUT/'inputs'/f'{method}.npz') as source:d={key:source[key] for key in source.files}
        gene_index={g:i for i,g in enumerate(d['genes'])};species_index={s:i for i,s in enumerate(d['species'])}
        chosen=set(range(min(5,len(rows))));chosen.update(rng.choice(len(rows),10,replace=False))
        chosen.update(i for i,r in enumerate(rows) if r['gene']=='SIRT6')
        for mode in ['interior','low_n','retained']:
            possible=[i for i,r in enumerate(rows) if (.1<float(r['lambda'])<.9 if mode=='interior' else int(r['n_species'])<=8 if mode=='low_n' else r['retained_OLS_FDR_hit']=='True')]
            chosen.update(possible[:3])
        for j in sorted(chosen):
            r=rows[j];i=gene_index[r['gene']];ix=[species_index[s] for s in r['species'].split(',')]
            points=[dict(species=str(d['scientific_names'][s]),x=float(d['x'][s]),y=float(d['y'][i,s]),
                         logMass=float(d['mass'][s]),lengthKb=float(d['length_kb'][i,s]) if method=='body' else 1.) for s in ix]
            cases.append(dict(method=method,gene=r['gene'],points=points,expected={k:float(r[k]) for k in ['MLS_beta','MLS_SE','p_value','lambda','GLS_R2','model_intercept']},
                              covariance=d['covariance'][np.ix_(ix,ix)].tolist()))
        # Reproduce previous OLS fits before interpreting PGLS attrition.
        oldfile=(ROOT/f"outputs/regulatory_windows_2026-10-22/tss_{method[5:]}/all_tested.tsv.gz" if method.startswith('peaks') else
                 ROOT/'outputs/gene_activity_tss5kb_2026-10-23/mass_corrected_all_genes.tsv' if method=='tss5kb' else
                 ROOT/f'outputs/activity_omics_comparison_2026-10-20/mass_adjusted/all_available/{method}_correlations.tsv')
        if oldfile.exists():
            opener=gzip.open if oldfile.suffix=='.gz' else open
            with opener(oldfile,'rt') as f:old={r['gene']:r for r in csv.DictReader(f,delimiter='\t')}
            assert old.keys()=={r['gene'] for r in rows},(method,'tested set changed')
            col='p_value' if method.startswith('peaks') else 'p'
            error=max(abs(float(r['OLS_p_value'])-float(old[r['gene']][col])) for r in rows)
            assert error<1e-7,(method,error)
            comparisons.append(dict(method=method,genes_compared=len(rows),maximum_OLS_p_difference=error,reference=str(oldfile)))
    (OUT/'validation_cases.json').write_text(json.dumps(cases))
    (OUT/'OLS_reproduction.json').write_text(json.dumps(comparisons,indent=2))
    print('Prepared',len(cases),'validation cases;',len(comparisons),'whole-scan OLS comparisons')

if __name__=='__main__':main()
