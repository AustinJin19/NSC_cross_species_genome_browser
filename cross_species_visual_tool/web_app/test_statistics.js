const assert=require('node:assert/strict');
const {linearFit}=require('./www/statistics.js');
const close=(a,b)=>assert.ok(Math.abs(a-b)<1e-12,`${a} != ${b}`);
let f=linearFit([{x:1,y:5},{x:2,y:7},{x:3,y:9}]);
close(f.slope,2);close(f.intercept,3);close(f.r,1);
f=linearFit([{x:1,y:9},{x:2,y:7},{x:3,y:5}]);close(f.r,-1);
f=linearFit([{x:1,y:1},{x:2,y:0},{x:3,y:1}]);close(f.r,0);
assert.equal(linearFit([{x:1,y:1},{x:2,y:2}]),null);
assert.equal(linearFit([{x:1,y:1},{x:1,y:2},{x:1,y:3}]),null);
assert.equal(linearFit([{x:1,y:2},{x:2,y:2},{x:3,y:2}]),null);
assert.equal(linearFit([{x:1,y:2},{x:2,y:NaN},{x:3,y:4}]),null);
console.log('Linear-fit checks passed');
const {pearsonP}=require('./www/statistics.js');
close(pearsonP(.47,4),.53);
close(pearsonP(-.47,4),.53);
close(pearsonP(.5,3),2/3);
close(pearsonP(0,10),1);
close(pearsonP(1,10),0);
assert.equal(pearsonP(.5,2),null);
console.log('Pearson p-value checks passed');
const {massAdjustedFit,logLengthMassAdjustedFit,logLengthAdjustedFit}=require('./www/statistics.js');
const fixtures=require('./test_mass_fit_fixtures.json');
for(const test of fixtures.cases){
 const fit=test.method==='massAdjustedFit'?massAdjustedFit(test.points):logLengthMassAdjustedFit(test.points);
 assert.ok(fit);
 for(const [key,expected] of Object.entries(test.expected)){
  assert.ok(Math.abs(fit[key]-expected)<1e-10,`${test.method} ${key}: ${fit[key]} != ${expected}`);
 }
}
const points=fixtures.cases[0].points;
assert.equal(massAdjustedFit(points.slice(0,3)),null);
assert.equal(logLengthMassAdjustedFit(points.slice(0,4)),null);
assert.equal(massAdjustedFit(points.map(p=>({...p,logMass:2}))),null);
assert.equal(massAdjustedFit(points.map(p=>({...p,logMass:p.x}))),null);
assert.equal(logLengthMassAdjustedFit(points.map(p=>({...p,logMass:Math.log10(p.lengthKb)}))),null);
assert.equal(logLengthMassAdjustedFit(points.map(p=>({...p,lengthKb:0}))),null);
assert.equal(massAdjustedFit(points.map(p=>({...p,logMass:NaN}))),null);
// An absent mass is irrelevant when the existing length-only option is used.
assert.ok(logLengthAdjustedFit(points.map(p=>({...p,logMass:NaN}))));
console.log('Mass-only and log-length + log-mass fits match independent SciPy/NumPy OLS on raw/log axes; invalid and collinear cases rejected.');
