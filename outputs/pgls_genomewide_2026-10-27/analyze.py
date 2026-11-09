"""PGLS estimated-lambda and matched OLS scans; one BH family per method."""
import argparse,csv,gzip,json,time
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import false_discovery_control
from pgls_batch import fit_batch

OUT=Path(__file__).resolve().parent
METHODS=['tss2kb','tss5kb','body','tss100kb','tss200kb','tss500kb','peaks20kb','peaks50kb','peaks100kb','peaks200kb','peaks500kb','peaks1000kb']
FIELDS=['method','gene','n_species','direction','MLS_beta','MLS_SE','p_value','BH_FDR','lambda','GLS_R2','df','model_intercept','mass_beta','length_beta','OLS_beta','OLS_p_value','OLS_BH_FDR','OLS_partial_r','retained_OLS_FDR_hit','species']

def write(path,rows,fields):
    path.parent.mkdir(exist_ok=True,parents=True)
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def scan(method):
    started=time.time()
    with np.load(OUT/'inputs'/f'{method}.npz') as source:d={key:source[key] for key in source.files}
    genes=d['genes'];y=d['y'];body=method=='body'
    n_genes,n_species=y.shape
    eligible=np.isfinite(y)
    if body:eligible&=np.isfinite(d['length_kb'])&(d['length_kb']>0)
    groups=defaultdict(list)
    for i,mask in enumerate(eligible):
        if mask.sum()>=5:groups[mask.tobytes()].append(i)
    rows=[];invalid=[]
    for counter,(mask_bytes,indexes) in enumerate(sorted(groups.items(),key=lambda kv:-len(kv[1]))):
        mask=np.frombuffer(mask_bytes,dtype=bool);ispecies=np.flatnonzero(mask);n=len(ispecies);C=d['covariance'][np.ix_(ispecies,ispecies)]
        # Limit transient QR arrays while sharing spectral decomposition across genes.
        for begin in range(0,len(indexes),2000):
            ids=np.asarray(indexes[begin:begin+2000]);Y=y[np.ix_(ids,ispecies)];g=len(ids)
            X=np.empty((g,n,4 if body else 3));X[:,:,0]=1;X[:,:,1]=d['x'][ispecies]
            if body:X[:,:,2]=np.log10(d['length_kb'][np.ix_(ids,ispecies)]);X[:,:,3]=d['mass'][ispecies]
            else:X[:,:,2]=d['mass'][ispecies]
            result=fit_batch(Y,X,C,'ML');ols=fit_batch(Y,X,C,0)
            for j,i in enumerate(ids):
                if not result['valid'][j] or not ols['valid'][j]:
                    invalid.append(dict(gene=str(genes[i]),n_species=n,PGLS_valid=bool(result['valid'][j]),OLS_valid=bool(ols['valid'][j])));continue
                beta=float(result['beta'][j,0]);obeta=float(ols['beta'][j,0]);ot=obeta/float(ols['se'][j,0]);df=int(result['df'][j])
                rows.append(dict(method=method,gene=str(genes[i]),n_species=n,direction='Positive' if beta>0 else 'Negative',
                    MLS_beta=beta,MLS_SE=float(result['se'][j,0]),p_value=float(result['p'][j]),**{'lambda':float(result['lambda_value'][j])},
                    GLS_R2=float(result['r_squared'][j]),df=df,model_intercept=float(result['intercept'][j]),
                    mass_beta=float(result['beta'][j,-1]),length_beta=float(result['beta'][j,1]) if body else None,
                    OLS_beta=obeta,OLS_p_value=float(ols['p'][j]),OLS_partial_r=float(ot/np.sqrt(ot*ot+df)),
                    species=','.join(str(v) for v in d['species'][ispecies])))
        if counter%100==0:print(method,'groups',counter+1,'/',len(groups),'fitted',len(rows),'elapsed',round(time.time()-started,1),flush=True)
    q=false_discovery_control([r['p_value'] for r in rows]);oq=false_discovery_control([r['OLS_p_value'] for r in rows])
    for r,a,b in zip(rows,q,oq):r['BH_FDR']=float(a);r['OLS_BH_FDR']=float(b);r['retained_OLS_FDR_hit']=bool(a<.05 and b<.05)
    rows.sort(key=lambda r:(r['p_value'],r['gene']))
    dest=OUT/'results'/method;write(dest/'all_tested.tsv.gz',rows,FIELDS)
    for cutoff,col in [('FDR05','BH_FDR'),('nominal_p05','p_value')]:
        hit=[r for r in rows if r[col]<.05];write(dest/f'{cutoff}.tsv',hit,FIELDS)
        for sign in ['Positive','Negative']:
            selected=[r for r in hit if r['direction']==sign]
            write(dest/f'{cutoff}_{sign.lower()}.tsv',selected,FIELDS)
            (dest/f'{cutoff}_{sign.lower()}_genes.txt').write_text(''.join(r['gene']+'\n' for r in selected))
    sig=[r for r in rows if r['BH_FDR']<.05];previous=[r for r in rows if r['OLS_BH_FDR']<.05]
    summary=dict(method=method,tested=len(rows),FDR_positive=sum(r['MLS_beta']>0 for r in sig),FDR_negative=sum(r['MLS_beta']<0 for r in sig),
        FDR_total=len(sig),OLS_FDR_total=len(previous),retained_OLS_FDR=sum(r['retained_OLS_FDR_hit'] for r in rows),
        new_PGLS_FDR=sum(r['BH_FDR']<.05 and r['OLS_BH_FDR']>=.05 for r in rows),lost_OLS_FDR=sum(r['BH_FDR']>=.05 and r['OLS_BH_FDR']<.05 for r in rows),
        nominal_positive=sum(r['p_value']<.05 and r['MLS_beta']>0 for r in rows),nominal_negative=sum(r['p_value']<.05 and r['MLS_beta']<0 for r in rows),
        minimum_FDR=min(r['BH_FDR'] for r in rows),minimum_n=min(r['n_species'] for r in rows),maximum_n=max(r['n_species'] for r in rows),
        lambda_zero=sum(r['lambda']<1e-6 for r in rows),lambda_one=sum(r['lambda']>1-1e-6 for r in rows),
        median_lambda=float(np.median([r['lambda'] for r in rows])),invalid_fits=len(invalid),elapsed_seconds=time.time()-started)
    (dest/'summary.json').write_text(json.dumps(summary,indent=2));(dest/'invalid_fits.json').write_text(json.dumps(invalid,indent=2))
    print('COMPLETE',json.dumps(summary),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('methods',nargs='+',choices=METHODS);args=p.parse_args()
    for m in args.methods:scan(m)
