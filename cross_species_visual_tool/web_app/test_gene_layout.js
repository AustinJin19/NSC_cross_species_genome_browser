const assert=require('node:assert/strict');
const {namedGene,packGenes}=require('./www/gene_layout.js');
for(const n of ['LOC100123','ENSHGLG00000012345','ENSG00000123456.2','XR_123.1','unknown',''])assert.equal(namedGene(n),false,n);
for(const n of ['APOE','APOC1','CEACAM20','Gm1234','C1orf54'])assert.equal(namedGene(n),true,n);
const models=[{gname:'APOE',start:100,end:200},{gname:'OTHER',start:110,end:210},{gname:'THIRD',start:120,end:220},{gname:'FOURTH',start:130,end:230},{gname:'APOC1',start:350,end:400},{gname:'LOC123',start:500,end:600}];
for(const position of [x=>x,x=>900-x]){
 const r=packGenes(models,'APOE',position,12,888,s=>s.length*5);
 assert.equal(new Set(r.placed.map(p=>p.lane)).size,4);
 assert.equal(r.filtered,1);assert.equal(r.rowCount,4);
 assert.equal(r.placed.length,5);
 assert.deepEqual(r.placed.map(p=>p.model.gname).sort(),['APOC1','APOE','FOURTH','OTHER','THIRD']);
 assert(r.placed.some(p=>p.model.gname==='APOE'));assert(r.placed.some(p=>p.model.gname==='APOC1'));
 for(const p of r.placed)assert(p.lane>=0&&p.lane<r.rowCount);
}
assert.equal(packGenes([{gname:'LOC123',start:1,end:10}],'LOC123',x=>x,0,900,s=>s.length*5).placed.length,1);
console.log('Adaptive gene layout checks passed');

function verifyLabels(r,left,right){
 for(let lane=0;lane<r.rowCount;lane++){
  const row=r.placed.filter(p=>p.lane===lane).sort((a,b)=>a.labelX-b.labelX);
  for(let i=0;i<row.length;i++){
   assert(row[i].labelX>=left-1e-8);
   assert(row[i].labelX+row[i].labelWidth<=right+1e-8);
   if(i)assert(row[i-1].labelX+row[i-1].labelWidth<row[i].labelX);
  }
 }
}
for(const n of [5,30,200])for(const flip of [false,true]){
 const models=Array.from({length:n},(_,i)=>({gname:'NAMED_GENE_'+i,start:450+i%3,end:480+i%7}));
 const r=packGenes(models,'NAMED_GENE_0',x=>flip?900-x:x,12,888,s=>s.length*5);
 assert.equal(r.placed.length,n);assert.equal(r.rowCount,n);verifyLabels(r,12,888);
 for(const p of r.placed){assert.equal(p.fontSize,8);assert.equal(p.labelWidth,p.model.gname.length*5);}
}
const disjoint=Array.from({length:8},(_,i)=>({gname:'GENE'+i,start:20+i*100,end:60+i*100}));
assert.equal(packGenes(disjoint,'GENE0',x=>x,12,888,s=>s.length*5).rowCount,1);
assert.equal(packGenes([],'GENE0',x=>x,12,888,s=>s.length*5).rowCount,0);
for(const p of packGenes(models,'APOE',x=>x,12,888,s=>s.length*5).placed)assert.equal(p.fontSize,8);
console.log('Label collision, compact packing and boundary checks passed');

const {intronArrows}=require('./www/gene_layout.js');
const transcript={strand:'+',ex_s:'100,200,',ex_e:'120,220,'};
for(const panelStrand of ['+','-']){
 const position=panelStrand==='+'?x=>x:x=>900-x;
 const arrows=intronArrows(transcript,panelStrand,position,0,900);
 assert(arrows.length>0);
 for(const a of arrows){
  assert.equal(a.direction,panelStrand==='+'?1:-1);
  assert(a.x>Math.min(position(121),position(200))&&a.x<Math.max(position(121),position(200)));
 }
}
assert.equal(intronArrows({strand:'+',ex_s:'100,',ex_e:'200,'},'+',x=>x,0,900).length,0);
assert.equal(intronArrows(transcript,'+',x=>x,300,900).length,0);
assert(intronArrows({...transcript,strand:'-'},'+',x=>x,0,900).every(a=>a.direction===-1));
console.log('Intron direction and clipping checks passed');
