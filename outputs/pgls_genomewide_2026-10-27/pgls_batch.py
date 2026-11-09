"""Batched Gaussian ML PGLS; same lambda model and t test as the browser.

Spectral whitening permits simultaneous fits with a shared taxon subset.
Each response has its own coefficients, variance and ML lambda. QR is used
instead of normal equations, including for gene-specific length covariates.
"""
import numpy as np
from scipy.stats import t


def fit_batch(Y, X, C, lambda_mode='ML'):
    """Y=(genes,species), X=(genes,species,parameters), including intercept."""
    Y=np.asarray(Y,float); X=np.asarray(X,float); C=np.asarray(C,float)
    g,n=Y.shape; k=X.shape[2]; df=n-k
    assert X.shape[:2]==Y.shape and df>0
    diag=np.diag(C); R=C/np.sqrt(diag[:,None]*diag[None,:])
    eig,U=np.linalg.eigh(R)
    if eig.min()<=1e-12:raise ValueError('Singular phylogenetic covariance')
    center=X[:,:,1:].mean(axis=1); scale=X[:,:,1:].std(axis=1)
    yc=Y.mean(axis=1); ys=Y.std(axis=1)
    valid=np.isfinite(X).all(axis=(1,2))&np.isfinite(Y).all(axis=1)&(scale>0).all(axis=1)&(ys>0)
    Z=X.copy(); Z[:,:,1:]=(X[:,:,1:]-center[:,None,:])/np.where(scale>0,scale,1)[:,None,:]
    y=(Y-yc[:,None])/np.where(ys>0,ys,1)[:,None]
    rotation=U.T/np.sqrt(diag)[None,:]
    xr=np.einsum('ij,gjk->gik',rotation,Z); yr=y@rotation.T
    log_d=np.log(diag).sum()

    def evaluate(ids,lam,details=False):
        ids=np.asarray(ids);lam=np.broadcast_to(np.asarray(lam,float),len(ids))
        eigen=1+lam[:,None]*(eig[None,:]-1)
        weights=1/np.sqrt(eigen)
        a=xr[ids]*weights[:,:,None];b=yr[ids]*weights
        q,r=np.linalg.qr(a,mode='reduced')
        rank_ok=(np.abs(np.diagonal(r,axis1=1,axis2=2))>1e-10*np.linalg.norm(a,axis=1)).all(axis=1)
        safe=r.copy();safe[~rank_ok]=np.eye(k)
        coef=np.linalg.solve(safe,np.einsum('gnk,gn->gk',q,b)[...,None])[...,0]
        err=b-np.einsum('gnk,gk->gn',a,coef);sse=np.einsum('gn,gn->g',err,err)
        okay=valid[ids]&rank_ok&np.isfinite(sse)&(sse>1e-24)
        likelihood=-.5*(n*(np.log(2*np.pi)+1+np.log(np.maximum(sse,1e-300)/n))+log_d+np.log(eigen).sum(axis=1))-n*np.log(np.maximum(ys[ids],1e-300))
        likelihood[~okay]=-np.inf
        if not details:return likelihood
        inverse=np.linalg.solve(safe,np.broadcast_to(np.eye(k),safe.shape))
        variance=(inverse**2).sum(axis=2)*(sse/df)[:,None]
        safe_scale=np.where(scale[ids]>0,scale[ids],1)
        betas=coef[:,1:]*ys[ids,None]/safe_scale
        se=np.sqrt(variance[:,1:])*ys[ids,None]/safe_scale
        intercept=yc[ids]+ys[ids]*coef[:,0]-(betas*center[ids]).sum(axis=1)
        null_intercept=(a[:,:,0]*b).sum(axis=1)/(a[:,:,0]**2).sum(axis=1)
        null_sse=((b-a[:,:,0]*null_intercept[:,None])**2).sum(axis=1)
        return dict(valid=okay,beta=betas,se=se,p=2*t.sf(np.abs(betas[:,0]/se[:,0]),df),
                    intercept=intercept,r_squared=1-sse/null_sse,log_likelihood=likelihood,
                    lambda_value=lam,df=np.full(len(ids),df))

    ids=np.arange(g)
    if lambda_mode!='ML':return evaluate(ids,float(lambda_mode),True)
    grid=np.linspace(0,1,21)
    scores=np.stack([evaluate(ids,v) for v in grid],axis=1)
    best_index=np.argmax(scores,axis=1);best_lambda=grid[best_index].copy();best_score=scores[ids,best_index].copy()
    left=np.column_stack([np.full(g,-np.inf),scores[:,:-1]])
    right=np.column_stack([scores[:,1:],np.full(g,-np.inf)])
    candidates=np.isfinite(scores)&(scores>=left)&(scores>=right)
    rows,columns=np.where(candidates)
    if len(rows):
        lo=np.maximum(0,(columns-1)/20);hi=np.minimum(1,(columns+1)/20)
        ratio=(np.sqrt(5)-1)/2
        c=hi-ratio*(hi-lo);d=lo+ratio*(hi-lo)
        fc=evaluate(rows,c);fd=evaluate(rows,d)
        for _ in range(40):
            use_left=fc>fd
            hi=np.where(use_left,d,hi);lo=np.where(use_left,lo,c)
            c=hi-ratio*(hi-lo);d=lo+ratio*(hi-lo)
            fc=evaluate(rows,c);fd=evaluate(rows,d)
        candidate_lambda=np.where(fc>fd,c,d);candidate_score=np.maximum(fc,fd)
        # Multiple local likelihood maxima for one gene are compared globally.
        for row,lam,score in zip(rows,candidate_lambda,candidate_score):
            if score>best_score[row]:best_score[row]=score;best_lambda[row]=lam
    return evaluate(ids,best_lambda,True)
