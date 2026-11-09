const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const app=path.resolve(__dirname,'../../cross_species_visual_tool/web_app');
const model=require(app+'/www/pgls.js'),tree=require(app+'/www/phylogeny.js');
const cases=JSON.parse(fs.readFileSync(__dirname+'/validation_cases.json'));
let maximum=0;
for(const c of cases){
 const result=model.fit(c.points,tree,{lambda:'ML',adjustLength:c.method==='body',adjustMass:true});
 assert.ok(result,c.method+':'+c.gene);
 const map={MLS_beta:'slope',MLS_SE:'slopeSE',p_value:'p',lambda:'lambda',GLS_R2:'rSquared',model_intercept:'modelIntercept'};
 for(const [expected,actual] of Object.entries(map)){
  const diff=Math.abs(c.expected[expected]-result[actual]);maximum=Math.max(maximum,diff);
  assert.ok(diff<2e-5,`${c.method}/${c.gene} ${actual}: ${diff}`);
 }
}
fs.writeFileSync(__dirname+'/browser_validation.json',JSON.stringify({fits_compared:cases.length,maximum_absolute_difference:maximum,tolerance:2e-5},null,2));
console.log(cases.length,'fits agree with production browser; maximum absolute difference',maximum);
