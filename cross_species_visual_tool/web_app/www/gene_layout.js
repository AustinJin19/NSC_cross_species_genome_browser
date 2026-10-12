'use strict';
// Pack complete gene models and labels into as few non-overlapping rows as possible.
function namedGene(symbol){
 const s=String(symbol||'').trim();
 return !!s && !/^(?:\.|-|NA|N\/A|unknown|unnamed)$/i.test(s)
  && !/^(?:LOC\d+|ENS[A-Z]*[GTP]\d+(?:\.\d+)?|(?:NM|NR|XM|XR)_\d+(?:\.\d+)?|gene[-_:].+)$/i.test(s);
}
function packGenes(models,gene,position,left,right,measure){
 const named=models.filter(m=>namedGene(m.gname)||m.gname.toUpperCase()===gene.toUpperCase());
 const candidates=named.map(m=>{
  const a=position(m.start),b=position(m.end),start=Math.max(left,Math.min(a,b)),end=Math.min(right,Math.max(a,b));
  const labelWidth=measure(m.gname),labelX=Math.max(left,Math.min(start,right-labelWidth));
  return {model:m,start:Math.min(start,labelX),end:Math.max(end,labelX+labelWidth),labelX,labelWidth,anchorX:start};
 }).sort((a,b)=>a.start-b.start||a.end-b.end||a.model.gname.localeCompare(b.model.gname));
 const laneEnds=[],placed=[];
 for(const c of candidates){
  // Left-to-right interval partitioning needs a new row only when all rows
  // overlap this model or its label. Keep labels at their measured size.
  let lane=laneEnds.findIndex(end=>end+6<=c.start);
  if(lane<0)lane=laneEnds.length;
  laneEnds[lane]=c.end;
  placed.push({...c,lane,fontSize:8});
 }
 return {placed,filtered:models.length-named.length,rowCount:laneEnds.length};
}
function intronArrows(model,panelStrand,position,left,right){
 const starts=model.ex_s.split(',').filter(Boolean).map(Number),ends=model.ex_e.split(',').filter(Boolean).map(Number);
 const exons=starts.map((a,i)=>[a,ends[i]]).filter(([a,b])=>Number.isFinite(a)&&Number.isFinite(b)&&b>a).sort((a,b)=>a[0]-b[0]);
 const arrows=[],direction=model.strand===panelStrand?1:-1;
 let previousEnd=exons[0]?.[1];
 for(let i=1;i<exons.length;i++){
  const [start,end]=exons[i];
  if(start>previousEnd){
   const a=position(previousEnd+1),b=position(start);
   const lo=Math.max(left,Math.min(a,b))+3,hi=Math.min(right,Math.max(a,b))-3;
   if(hi>=lo){
    const count=Math.max(1,Math.floor((hi-lo)/30)+1);
    for(let j=0;j<count;j++)arrows.push({x:lo+(j+.5)*(hi-lo)/count,direction});
   }
  }
  previousEnd=Math.max(previousEnd,end);
 }
 return arrows;
}
if(typeof module!=='undefined')module.exports={namedGene,packGenes,intronArrows};
