/* PGLS with Brownian shared-branch covariance and optional ML Pagel's lambda.
 * Cholesky whitening + reorthogonalized QR; no diagonal jitter or OLS fallback.
 * Species medians are the observations. The response and covariates are chosen
 * by the calling panel. Coefficient t tests condition on the fitted lambda.
 */
const NSCPGLS = (() => {
 const statistics = typeof NSCStatistics!=='undefined'?NSCStatistics:require('./statistics.js');
 const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
 const mean=a=>a.reduce((s,v)=>s+v,0)/a.length;
 const name=p=>p.sp?.latin||p.species;
 function covarianceFor(points,tree){
  if(!tree||!Array.isArray(tree.species)||!Array.isArray(tree.covariance))return null;
  const names=points.map(name);
  if(new Set(names).size!==names.length||new Set(tree.species).size!==tree.species.length)return null;
  const ix=names.map(s=>tree.species.indexOf(s));
  if(ix.some(i=>i<0))return null;
  const C=ix.map(i=>ix.map(j=>tree.covariance[i]?.[j]));
  if(C.some((row,i)=>row.some((v,j)=>!Number.isFinite(v)||v<0||Math.abs(v-C[j][i])>1e-9)||row[i]<=0))return null;
  return C;
 }
 function cholesky(C){
  const n=C.length,L=Array.from({length:n},()=>Array(n).fill(0));
  for(let i=0;i<n;i++)for(let j=0;j<=i;j++){
   let v=C[i][j];for(let k=0;k<j;k++)v-=L[i][k]*L[j][k];
   if(i===j){if(v<=1e-12||!Number.isFinite(v))return null;L[i][j]=Math.sqrt(v);}
   else L[i][j]=v/L[j][j];
  }
  return L;
 }
 function forward(L,b){
  const x=[];for(let i=0;i<b.length;i++){let v=b[i];for(let j=0;j<i;j++)v-=L[i][j]*x[j];x.push(v/L[i][i]);}return x;
 }
 function backward(R,b){
  const x=Array(b.length).fill(0);for(let i=b.length-1;i>=0;i--){let v=b[i];for(let j=i+1;j<b.length;j++)v-=R[i][j]*x[j];x[i]=v/R[i][i];}return x;
 }
 function qrSolve(columns,y){
  const k=columns.length,Q=[],R=Array.from({length:k},()=>Array(k).fill(0));
  for(let j=0;j<k;j++){
   let v=[...columns[j]],original=Math.sqrt(dot(v,v));
   for(let pass=0;pass<2;pass++)for(let i=0;i<j;i++){
    const z=dot(Q[i],v);R[i][j]+=z;v=v.map((x,h)=>x-z*Q[i][h]);
   }
   const norm=Math.sqrt(dot(v,v));
   if(norm<=1e-10*original||!Number.isFinite(norm))return null;
   R[j][j]=norm;Q.push(v.map(x=>x/norm));
  }
  const beta=backward(R,Q.map(q=>dot(q,y)));
  const residual=y.map((v,i)=>v-columns.reduce((s,c,j)=>s+c[i]*beta[j],0));
  const inverse=Array.from({length:k},(_,j)=>backward(R,Array.from({length:k},(_,i)=>+(i===j))));
  return {beta,sse:dot(residual,residual),varianceDiagonal:beta.map((_,i)=>inverse.reduce((s,col)=>s+col[i]**2,0))};
 }
 function fit(points,tree,{lambda='ML',adjustLength=false,adjustMass=false}={}){
  const n=points.length,k=2+Number(adjustLength)+Number(adjustMass),df=n-k;
  if(n<3||df<1||points.some(p=>![p.x,p.y].every(Number.isFinite)))return null;
  if(lambda!=='ML'&&(!Number.isFinite(lambda)||lambda<0||lambda>1))return null;
  if(adjustLength&&points.some(p=>!Number.isFinite(p.lengthKb)||p.lengthKb<=0))return null;
  if(adjustMass&&points.some(p=>!Number.isFinite(p.logMass)))return null;
  let C=covarianceFor(points,tree);if(!C)return null;
  const cScale=mean(C.map((r,i)=>r[i]));C=C.map(r=>r.map(v=>v/cScale));
  const rawColumns=[points.map(p=>p.x)];
  if(adjustLength)rawColumns.push(points.map(p=>Math.log10(p.lengthKb)));
  if(adjustMass)rawColumns.push(points.map(p=>p.logMass));
  const centers=rawColumns.map(mean),scales=rawColumns.map((c,j)=>Math.sqrt(mean(c.map(v=>(v-centers[j])**2))));
  if(scales.some(s=>s<=0||!Number.isFinite(s)))return null;
  const X=[Array(n).fill(1),...rawColumns.map((c,j)=>c.map(v=>(v-centers[j])/scales[j]))];
  const my=mean(points.map(p=>p.y)),sy=Math.sqrt(mean(points.map(p=>(p.y-my)**2)));
  if(sy<=0||!Number.isFinite(sy))return null;
  const y=points.map(p=>(p.y-my)/sy);
  function evaluate(lam){
   const L=cholesky(C.map((r,i)=>r.map((v,j)=>i===j?v:lam*v)));if(!L)return null;
   const columns=X.map(col=>forward(L,col)),response=forward(L,y),sol=qrSolve(columns,response);
   if(!sol||sol.sse<=1e-24||!Number.isFinite(sol.sse))return null;
   const logDet=2*L.reduce((s,row,i)=>s+Math.log(row[i]),0);
   const logLikelihood=-.5*(n*(Math.log(2*Math.PI)+1+Math.log(sol.sse/n))+logDet)-n*Math.log(sy);
   const nullIntercept=dot(columns[0],response)/dot(columns[0],columns[0]);
   const nullSSE=response.reduce((s,v,i)=>s+(v-nullIntercept*columns[0][i])**2,0);
   return {...sol,lambda:lam,logLikelihood,rSquared:1-sol.sse/nullSSE};
  }
  let best=null;
  const consider=v=>{if(v&&(!best||v.logLikelihood>best.logLikelihood))best=v;return v;};
  if(lambda==='ML'){
   // Scan the entire interval and refine each local maximum, including boundary
   // neighborhoods. Evaluate exact endpoints so zero signal is a valid result.
   const grid=Array.from({length:21},(_,i)=>consider(evaluate(i/20)));
   const score=v=>v?v.logLikelihood:-Infinity;
   for(let i=0;i<grid.length;i++){
    if(score(grid[i])<score(grid[i-1])||score(grid[i])<score(grid[i+1]))continue;
    let a=Math.max(0,(i-1)/20),b=Math.min(1,(i+1)/20),ratio=(Math.sqrt(5)-1)/2;
    let c=b-ratio*(b-a),d=a+ratio*(b-a),fc=consider(evaluate(c)),fd=consider(evaluate(d));
    for(let iter=0;iter<80&&b-a>1e-8;iter++){
     if(score(fc)>score(fd)){b=d;d=c;fd=fc;c=b-ratio*(b-a);fc=consider(evaluate(c));}
     else {a=c;c=d;fc=fd;d=a+ratio*(b-a);fd=consider(evaluate(d));}
    }
   }
  }else best=evaluate(lambda);
  if(!best)return null;
  const coefficients=best.beta.slice(1).map((b,j)=>b*sy/scales[j]);
  const modelIntercept=my+sy*best.beta[0]-coefficients.reduce((s,b,j)=>s+b*centers[j],0);
  const slope=coefficients[0],slopeSE=Math.sqrt(best.sse/df*best.varianceDiagonal[1])*sy/scales[0];
  const t=slope/slopeSE,rForP=t/Math.hypot(t,Math.sqrt(df));
  let index=1;
  const lengthSlope=adjustLength?coefficients[index++]:undefined;
  const massSlope=adjustMass?coefficients[index++]:undefined;
  return {slope,slopeSE,t,p:statistics.pearsonP(rForP,df+2),df,n,
   modelIntercept,intercept:my+sy*best.beta[0]-slope*centers[0],lengthSlope,massSlope,
   lambda:best.lambda,lambdaEstimated:lambda==='ML',lambdaAtBoundary:best.lambda<1e-6||best.lambda>1-1e-6,
   rSquared:best.rSquared,logLikelihood:best.logLikelihood,
   species:points.map(name)};
 }
 return {fit,covarianceFor};
})();
if(typeof module!=='undefined')module.exports=NSCPGLS;
