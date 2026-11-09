const assert=require('node:assert/strict');
const pgls=require('./www/pgls.js');
const tree=require('./www/phylogeny.js');
const stats=require('./www/statistics.js');
const fixtures=require('./test_pgls_fixtures.json');
const close=(a,b,tol=1e-8)=>assert.ok(Number.isFinite(a)&&Math.abs(a-b)<tol,`${a} != ${b}`);
assert.deepEqual(tree,require('./data/phylogenetic_covariance.json'));
for(const test of fixtures.cases){
 const f=pgls.fit(test.points,tree,test.options);assert.ok(f,test.name);
 const tol=test.options.lambda==='ML'?2e-5:1e-8;
 for(const [key,value] of Object.entries(test.expected)){
  try{close(f[key],value,tol);}catch(e){throw Error(`${test.name} ${key}: ${e.message}`);}
 }
 const permuted=pgls.fit([...test.points].reverse(),tree,test.options);
 close(permuted.slope,f.slope,1e-7);close(permuted.p,f.p,1e-7);close(permuted.lambda,f.lambda,1e-6);
 console.log(`${test.name}: lambda=${f.lambda.toFixed(5)}, beta=${f.slope.toFixed(5)}, p=${f.p.toFixed(6)} matches R caper.`);
}
const points=fixtures.cases[0].points;
const ordinary=stats.linearFit(points),zero=pgls.fit(points,tree,{lambda:0});
close(zero.slope,ordinary.slope);close(zero.intercept,ordinary.intercept);close(zero.p,ordinary.p);close(zero.rSquared,ordinary.r**2);
const adjusted=stats.logLengthMassAdjustedFit(points),adjustedZero=pgls.fit(points,tree,{lambda:0,adjustLength:true,adjustMass:true});
for(const key of ['slope','intercept','modelIntercept','lengthSlope','massSlope','p'])close(adjustedZero[key],adjusted[key]);
const scaled={...tree,covariance:tree.covariance.map(row=>row.map(v=>v*1000))};
const a=pgls.fit(points,tree),b=pgls.fit(points,scaled);close(a.slope,b.slope,1e-7);close(a.p,b.p,1e-7);close(a.lambda,b.lambda,1e-6);
assert.equal(pgls.fit(points.slice(0,2),tree),null);
assert.equal(pgls.fit(points.slice(0,4),tree,{adjustLength:true,adjustMass:true}),null);
assert.equal(pgls.fit(points.map(p=>({...p,x:1})),tree),null);
assert.equal(pgls.fit(points.map(p=>({...p,y:1})),tree),null);
assert.equal(pgls.fit(points.map(p=>({...p,logMass:p.x})),tree,{adjustMass:true}),null);
assert.equal(pgls.fit(points.map(p=>({...p,lengthKb:0})),tree,{adjustLength:true}),null);
assert.equal(pgls.fit(points.map(p=>({...p,logMass:NaN})),tree,{adjustMass:true}),null);
assert.equal(pgls.fit([{...points[0],species:'Unknown'},...points.slice(1)],tree),null);
assert.equal(pgls.fit([points[0],...points.slice(0,-1)],tree),null);
assert.equal(pgls.fit(points,tree,{lambda:2}),null);
assert.equal(pgls.fit(points,{...tree,covariance:tree.covariance.map(r=>r.map(()=>1))},{lambda:1}),null);
console.log('PGLS matches OLS at lambda=0; invariant to species ordering and time units; rejects missing taxa, singular covariance and invalid/collinear data.');
