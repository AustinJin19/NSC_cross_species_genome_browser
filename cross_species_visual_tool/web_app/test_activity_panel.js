// Exercise both real panel renderers with synthetic API data and a minimal DOM.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
class Element {
 constructor(){this.children=[];this.options=[{},{}];this.value='';this.checked=false;this.textContent='';}
 append(...nodes){this.children.push(...nodes);}
 replaceChildren(...nodes){this.children=nodes;}
 setAttribute(){}
}
const statusIndex = () => 5;
const elements = new Map(), get = id => {
 if(!elements.has(id))elements.set(id,new Element());return elements.get(id);
};
for(const prefix of ['activity','strength']){
 get(prefix+'Scope').value='all';get(prefix+'Axis').value='years';get(prefix+'Scale').value='raw';
 get(prefix+'Strategy').value=prefix==='activity'?'body':'peaks100kb';
}
get('strengthMeasure').value='count';
const catalog=['short','half','full','multi'].map((id,i)=>({id,name:id,latin:id,color:'#333'}));
const record=(id,median,fraction,locus='G|chr1:3000-4000:+')=>({species:id,median,locus,
 samples:{sample:median},unit:'CPM_signal_bp',available_fraction:fraction,region_eligible:fraction>=.5,positive_signal_fraction:.5,signal_eligible:true});
const records=[record('short',100,.49),record('half',20,.5),record('full',30,1),
 record('multi',200,.4),record('multi',10,1,'G|chr1:5000-6000:+')];
const fits=[];
const correctedFits=[];
const pglsFits=[];
const lifespans=Object.fromEntries(catalog.map((s,i)=>[s.id,{years:5+i,log10_years:Math.log10(5+i),body_mass_g:s.id==='full'?null:100}]));
const context={console,document:{getElementById:get,createElement:()=>new Element(),createElementNS:()=>new Element()},
 state:{gene:'G',catalog,species:catalog.map(s=>s.id)},window:{},
 NSCPhylogeny:{species:catalog.map(s=>s.latin)},
 NSCPGLS:{fit:(points,tree,options)=>{pglsFits.push({points,options});return null;}},
 NSCStatistics:{linearFit:points=>{fits.push(points);return null;},logLengthAdjustedFit:points=>{fits.push(points);return null;},
  massAdjustedFit:points=>{correctedFits.push({method:'mass',points});return null;},
  logLengthMassAdjustedFit:points=>{correctedFits.push({method:'length+mass',points});return null;}},
 fetch:async url=>({ok:true,json:async()=>({gene:'G',strategy:url.split('strategy=')[1],records,
 lifespans})})};
vm.runInNewContext(fs.readFileSync(__dirname+'/www/activity.js','utf8'),context);
setImmediate(async()=>{
 assert.equal(fits.length,2);
 for(const points of fits)assert.deepEqual(Array.from(points,p=>p.sp.id),['half','full']);
 for(const prefix of ['activity','strength']){
  const rows=get(prefix+'Rows').children.map(tr=>tr.children.map(td=>td.textContent));
  assert.match(rows[0][statusIndex(prefix)],/Excluded: <50%/);
  assert.equal(rows[1][statusIndex(prefix)],'Included');
  assert.equal(rows[1][4],'50%');
  assert.match(rows[3][statusIndex(prefix)],/Excluded: <50%.*highest-activity locus/);
  assert.match(get(prefix+'Status').textContent,/2 species included/);
  get(prefix+'Scale').value='log10';get(prefix+'Axis').value='log10_years';
  get(prefix+'Scale').onchange();
 }
 for(const points of fits.slice(2)){
  assert.deepEqual(Array.from(points,p=>p.sp.id),['half','full']);
  assert.equal(points[0].y,Math.log10(20));
 }
 for(const [prefix,method] of [['activity','length+mass'],['strength','mass']]){
  assert.equal(get(prefix+'MassCorrection').checked,false);
  get(prefix+'MassCorrection').checked=true;get(prefix+'MassCorrection').onchange();
  const fit=correctedFits.at(-1);
  assert.equal(fit.method,method);
  assert.deepEqual(Array.from(fit.points,p=>p.sp.id),['half']);
  assert.equal(fit.points[0].logMass,2);
  assert.equal(get(prefix+'Rows').children[2].children[statusIndex(prefix)].textContent,'Body mass missing or nonpositive');
  assert.match(get(prefix+'Status').textContent,/log₁₀\(body mass in g\)/);
  get(prefix+'MassCorrection').checked=false;get(prefix+'MassCorrection').onchange();
  assert.deepEqual(Array.from(fits.at(-1),p=>p.sp.id),['half','full']);
 }
 for(const mass of [0,-10]){
  lifespans.full.body_mass_g=mass;
  get('strengthMassCorrection').checked=true;get('strengthMassCorrection').onchange();
  assert.deepEqual(Array.from(correctedFits.at(-1).points,p=>p.sp.id),['half']);
 }
 // Window activity and CPM-in-peaks use mass alone; gene body keeps log length.
 get('activityMassCorrection').checked=true;get('activityStrategy').value='tss100kb';
 await get('activityStrategy').onchange();
 assert.equal(correctedFits.at(-1).method,'mass');
 get('strengthMeasure').value='signal';await get('strengthMeasure').onchange();
 assert.equal(correctedFits.at(-1).method,'mass');
 assert.match(get('strengthStatus').textContent,/Signal is integrated only over the union/);
 // PGLS is independently selectable in each panel and receives the same
 // scale, coverage, highest-locus and covariate filters as ordinary models.
 lifespans.full.body_mass_g=150;
 get('activityStrategy').value='body';await get('activityStrategy').onchange();
 for(const prefix of ['activity','strength']){
  get(prefix+'MassCorrection').checked=true;
  get(prefix+'Regression').value='pgls_ml';get(prefix+'Regression').onchange();
  let f=pglsFits.at(-1);
  assert.equal(f.options.lambda,'ML');assert.equal(f.options.adjustLength,prefix==='activity');
  assert.equal(f.options.adjustMass,true);
  assert.deepEqual(Array.from(f.points,p=>p.sp.id),['half','full']);
  assert.equal(f.points[0].y,Math.log10(20));
  assert.match(get(prefix+'Status').textContent,/PGLS uses shared branch lengths/);
  assert.match(get(prefix+'Chart').children.at(-1).textContent,/PGLS requires/);
  assert.doesNotMatch(get(prefix+'Chart').children.at(-1).textContent,/Pearson/);
  get(prefix+'Regression').value='pgls_bm';get(prefix+'Regression').onchange();
  assert.equal(pglsFits.at(-1).options.lambda,1);
  // A missing tree tip must be explained in the table, never silently fitted.
  context.NSCPhylogeny.species=['half'];get(prefix+'Regression').onchange();
  assert.deepEqual(Array.from(pglsFits.at(-1).points,p=>p.sp.id),['half']);
  assert.equal(get(prefix+'Rows').children[2].children[statusIndex(prefix)].textContent,'Species absent from TimeTree');
  get(prefix+'Regression').value='ols';get(prefix+'Regression').onchange();
  assert.deepEqual(Array.from(correctedFits.at(-1).points,p=>p.sp.id),['half','full']);
  context.NSCPhylogeny.species=catalog.map(s=>s.latin);
 }
 // Human stays in the tables but is excluded from every fit; macaque is eligible.
 for(const [id,latin,years,mass] of [['Human','Homo sapiens',122.5,62040],['Macaque','Macaca mulatta',40,8235]]){
  catalog.push({id,name:id,latin,color:'#333'});
  records.push(record(id,50,1));
  lifespans[latin]={years,log10_years:Math.log10(years),body_mass_g:mass};
 }
 context.NSCPhylogeny.species=catalog.map(s=>s.latin);
 for(const [prefix,strategy,measure] of [['activity','body',null],['activity','tss2kb',null],['strength','peaks200kb','count'],['strength','peaks200kb','signal']]){
  get(prefix+'Strategy').value=strategy;
  if(measure)get('strengthMeasure').value=measure;
  await get(prefix+'Strategy').onchange();
  for(const scale of ['raw','log10'])for(const regression of ['ols','pgls_ml','pgls_bm'])for(const mass of [false,true]){
   get(prefix+'Scale').value=scale;get(prefix+'Regression').value=regression;
   get(prefix+'MassCorrection').checked=mass;get(prefix+'Scale').onchange();
   const points=regression==='ols'?(mass?correctedFits.at(-1).points:fits.at(-1)):pglsFits.at(-1).points;
   for(const id of ['Human','Macaque']){
    assert.equal(points.some(p=>p.sp.id===id),id==='Macaque',`${id}: ${prefix} ${measure} ${scale} ${regression} mass=${mass}`);
    const row=get(prefix+'Rows').children.find(tr=>tr.children[0].textContent===id);
    assert.equal(row.children[statusIndex(prefix)].textContent,id==='Human'?'Excluded from fit':'Included');
   }
  }
 }
 // Historical positive-coverage fields must not filter any live gene-activity fit.
 records[1].positive_signal_fraction=.4999;records[1].signal_eligible=false;
 records[2].positive_signal_fraction=.5;records[2].signal_eligible=true;
 records[3].available_fraction=1;records[3].region_eligible=true;
 records[3].positive_signal_fraction=.1;records[3].signal_eligible=false;
 for(const strategy of ['body','tss2kb','tss5kb','tss100kb','tss200kb','tss500kb']){
  get('activityStrategy').value=strategy;await get('activityStrategy').onchange();
  for(const regression of ['ols','pgls_ml','pgls_bm'])for(const scale of ['raw','log10']){
   get('activityRegression').value=regression;get('activityScale').value=scale;
   get('activityMassCorrection').checked=false;get('activityScale').onchange();
   const points=regression==='ols'?fits.at(-1):pglsFits.at(-1).points;
   assert.deepEqual(Array.from(points,p=>p.sp.id),['half','full','multi','Macaque']);
   assert.equal(get('activityRows').children[1].children[5].textContent,'Included');
   assert.equal(get('activityRows').children[2].children[5].textContent,'Included');
   assert.match(get('activityRows').children[3].children[5].textContent,/Included.*highest-activity locus/);
   assert.equal(points.find(p=>p.sp.id==='multi').y,scale==='raw'?200:Math.log10(200));
   assert.equal(get('activityRows').children[0].children.length,6);
   assert.equal(get('activitySamples').children[0].children.length,8);
   assert.doesNotMatch(get('activityStatus').textContent,/positive signal|CPM > 0/);
  }
 }
 // Unknown or absent positive coverage is also irrelevant to eligibility.
 records[2].positive_signal_fraction=null;records[2].signal_eligible=false;
 get('activityScale').onchange();
 assert.equal(get('activityRows').children[2].children[5].textContent,'Included');
 assert.deepEqual(Array.from(pglsFits.at(-1).points,p=>p.sp.id),['half','full','multi','Macaque']);
 for(const record of records){delete record.positive_signal_fraction;delete record.signal_eligible;}
 get('activityScale').onchange();
 assert.deepEqual(Array.from(pglsFits.at(-1).points,p=>p.sp.id),['half','full','multi','Macaque']);
 // Zero total activity remains valid on the raw scale, but log10(0) is undefined.
 records[2].median=0;
 for(const regression of ['ols','pgls_ml','pgls_bm'])for(const scale of ['raw','log10']){
  get('activityRegression').value=regression;get('activityScale').value=scale;get('activityScale').onchange();
  const points=regression==='ols'?fits.at(-1):pglsFits.at(-1).points;
  assert.equal(points.some(p=>p.sp.id==='full'),scale==='raw');
  assert.equal(get('activityRows').children[2].children[5].textContent,scale==='raw'?'Included':'Nonpositive activity (log undefined)');
 }
 records[2].median=30;
 // Peak-count / CPM-in-peaks fits retain their existing behavior.
 get('strengthMassCorrection').checked=false;get('strengthRegression').value='ols';
 for(const measure of ['count','signal']){
  get('strengthMeasure').value=measure;await get('strengthMeasure').onchange();
  assert.deepEqual(Array.from(fits.at(-1),p=>p.sp.id),['half','full','multi','Macaque']);
 }
 console.log('Sparse, unknown and absent positive-signal coverage no longer exclude species in any gene-activity window or regression; raw zeros remain valid and log zeros are excluded.');
 console.log('Human is excluded and macaque included in gene-body, TSS, peak-count and CPM-in-peaks fits on raw/log scales, with optional mass correction and OLS/PGLS.');
 console.log('Both panels exclude <50%, retain exactly 50%, preserve highest-locus selection, and fit only eligible species on raw/log scales.');
 console.log('Independent mass toggles select the correct models, exclude invalid mass only when enabled, and cover gene-body, TSS, peak-count and CPM-in-peaks methods.');
 console.log('PGLS selectors preserve panel filters and covariates, distinguish fixed/estimated lambda, and report missing tree species without an OLS fallback.');
});
