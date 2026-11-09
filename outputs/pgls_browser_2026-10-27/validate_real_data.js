// Run the production panel renderer with saved local API responses, then
// compare the exact filtered points and coefficients independently in R.
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const app=path.resolve(__dirname,'../../cross_species_visual_tool/web_app');
const tree=require(app+'/www/phylogeny.js'),pgls=require(app+'/www/pgls.js');
const statistics=require(app+'/www/statistics.js');
const catalog=JSON.parse(fs.readFileSync(__dirname+'/catalog.json'));
class Element {constructor(){this.children=[];this.options=[{},{}];this.value='';this.checked=false;this.textContent='';}append(...x){this.children.push(...x);}replaceChildren(...x){this.children=x;}setAttribute(){}}
const elements=new Map(), get=id=>{if(!elements.has(id))elements.set(id,new Element());return elements.get(id);};
for(const p of ['activity','strength']){
 get(p+'Scope').value='all';get(p+'Axis').value='log10_years';get(p+'Scale').value='log10';
 get(p+'Strategy').value=p==='activity'?'tss5kb':'peaks50kb';
 get(p+'Regression').value='pgls_ml';get(p+'MassCorrection').checked=true;get(p+'ExcludeBowhead').checked=true;
}
get('strengthMeasure').value='count';
const fits=[];let name=null;
const context={console,document:{getElementById:get,createElement:()=>new Element(),createElementNS:()=>new Element()},
 state:{gene:'SIRT6',catalog,species:catalog.map(s=>s.id)},window:{},NSCPhylogeny:tree,NSCStatistics:statistics,
 NSCPGLS:{fit:(points,t,options)=>{const result=pgls.fit(points,t,options);if(name)fits.push({name,options,points:points.map(p=>({species:p.sp.latin,x:p.x,y:p.y,lengthKb:p.lengthKb,logMass:p.logMass})),result});return result;}},
 fetch:async url=>({ok:true,json:async()=>JSON.parse(fs.readFileSync(__dirname+'/SIRT6_'+url.split('strategy=')[1]+'.json'))})};
vm.runInNewContext(fs.readFileSync(app+'/www/activity.js','utf8'),context);
setImmediate(async()=>{
 fits.length=0;
 for(const strategy of ['tss5kb','peaks50kb','body']){
 name=null;const prefix=strategy==='peaks50kb'?'strength':'activity';get(prefix+'Strategy').value=strategy;await get(prefix+'Strategy').onchange();
 for(const mode of ['pgls_ml','pgls_bm'])for(const mass of [false,true]){
  name=`SIRT6_${strategy}_${mode}_mass_${mass}`;get(prefix+'Regression').value=mode;get(prefix+'MassCorrection').checked=mass;get(prefix+'Regression').onchange();
 }
 }
 const unique=fits;
 fs.writeFileSync(__dirname+'/real_data_fits.json',JSON.stringify(unique,null,2)+'\n');
 for(const f of unique)console.log(f.name,JSON.stringify(f.result));
});
