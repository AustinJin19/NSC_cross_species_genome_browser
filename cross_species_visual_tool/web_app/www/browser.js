'use strict';
const $=id=>document.getElementById(id);
const state={catalog:[],species:[],panels:[],nextPanel:1,bins:700,gene:'HAS2',window:5000,pan:0,data:null,request:0,loading:false};
const fmt=n=>Math.round(n).toLocaleString('en-US');
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const svgNS='http://www.w3.org/2000/svg';
const tip=document.createElement('div');tip.id='tooltip';tip.hidden=true;document.body.append(tip);
function svgEl(tag,attrs={},text){const e=document.createElementNS(svgNS,tag);for(const[k,v]of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;return e;}
function add(svg,tag,attrs,text){const e=svgEl(tag,attrs,text);svg.append(e);return e;}
function sizeLabel(bp){return Math.abs(bp)>=1000?`${+(bp/1000).toFixed(1)} kb`:`${Math.round(bp)} bp`;}
async function json(url,options){const r=await fetch(url,options);const d=await r.json();if(!r.ok)throw Error(d.error||'Request failed');return d;}
function orderPanels(panels){const order=new Map(state.catalog.map((s,i)=>[s.id,i]));return panels.sort((a,b)=>order.get(a.species||a.id)-order.get(b.species||b.id));}
function panelSpec(){return orderPanels(state.panels).filter(p=>state.species.includes(p.species)).map(p=>({key:p.key,species:p.species,sample:p.sample}));}
function annotationLabel(a){return `${a.gene}: ${a.samples}/${a.totalSamples} samples across ${a.species}/${a.totalSpecies} species have this gene annotation`;}
function newPanel(meta,sample){return {key:'panel-'+state.nextPanel++,species:meta.id,sample:sample||meta.default,min:null,max:null};}
async function load(){
 window.NSCOmics?.show(state.gene);
 const request=++state.request;state.loading=true;$('tracks').setAttribute('aria-busy','true');tip.hidden=true;$('exportPdf').disabled=true;$('status').textContent='Reading local genomic intervals…';
 try{const data=await json('/api/browse',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({gene:state.gene,window:state.window,pan:state.pan,bins:state.bins,panels:panelSpec()})});if(request!==state.request)return;state.data=data;state.loading=false;$('tracks').setAttribute('aria-busy','false');render();}
 catch(e){if(request!==state.request)return;state.loading=false;$('tracks').setAttribute('aria-busy','false');state.data=null;$('tracks').replaceChildren();$('ruler').replaceChildren();$('status').textContent=e.message;$('summary').textContent='Query failed';}
}
function render(){
 const d=state.data;if(!d)return;orderPanels(d.panels);
 $('locusTitle').textContent=d.gene;$('summary').textContent=`${new Set(d.panels.filter(p=>!p.error).map(p=>p.id)).size} species · ${d.panels.length} panels · ${fmt(d.window*2)} bp · ${Math.min(d.bins,d.window*2+1)} bins`;
 if($('gene').value.trim().toUpperCase()===d.gene)$('annotationCount').textContent=annotationLabel(d.annotation);
 $('speciesCount').textContent=state.species.length;
 const good=d.panels.filter(p=>!p.error), finite=good.flatMap(p=>p.signal.map(v=>v.value).filter(Number.isFinite));
 const globalMax=finite.reduce((a,b)=>Math.max(a,b),.001), globalMin=finite.reduce((a,b)=>Math.min(a,b),0);
 $('exportPdf').disabled=state.loading||!finite.length;
 $('status').textContent=state.loading?'Reading local genomic intervals…':!d.panels.length?'Select at least one species to explore.':good.some(p=>p.ambiguous)?'Multiple loci share this symbol in one or more species; the longest locus is shown.':'';
 $('tracks').replaceChildren();
 const width=900,pad=12,x=rel=>pad+(rel+d.window)/(2*d.window)*(width-2*pad);
 const ruler=svgEl('svg',{viewBox:`0 0 ${width} 37`,preserveAspectRatio:'none',width:'100%',height:37});
 for(let i=0;i<=4;i++){const rel=-d.window+i*d.window/2;add(ruler,'text',{x:x(rel),y:22,fill:'#6c757d','font-size':9,'font-family':'monospace','text-anchor':i===0?'start':i===4?'end':'middle'},sizeLabel(rel+d.pan));add(ruler,'line',{x1:x(rel),x2:x(rel),y1:29,y2:37,stroke:'#dee2e6'});}
 $('ruler').replaceChildren(ruler);
 for(const p of d.panels){
  const meta=state.catalog.find(s=>s.id===p.id), row=document.createElement('div');row.className='track';row.dataset.panelKey=p.panelKey;const config=state.panels.find(c=>c.key===p.panelKey);if(!config)continue;
  const label=document.createElement('div');label.className='track-label';label.innerHTML=`<div class="track-name"><span class="dot" style="background:${meta.color}"></span>${esc(meta.name)}</div><em>${esc(meta.latin)}</em>`;
  const select=document.createElement('select');select.setAttribute('aria-label',`${meta.name} replicate ${p.panelKey}`);for(const sample of meta.samples){const o=new Option(sample,sample);o.selected=sample===config.sample;select.add(o);}select.onchange=()=>{config.sample=select.value;load();};label.append(select);
  const actions=document.createElement('div');actions.className='panel-actions';
  const duplicate=document.createElement('button');duplicate.type='button';duplicate.textContent='+ Add panel';duplicate.setAttribute('aria-label',`Add ${meta.name} panel`);duplicate.onclick=()=>{if(state.panels.length>=40){$('status').textContent='Up to 40 panels are supported.';return;}const used=state.panels.filter(c=>c.species===p.id).map(c=>c.sample),sample=meta.samples.find(s=>!used.includes(s))||config.sample;state.panels.splice(state.panels.indexOf(config)+1,0,newPanel(meta,sample));load();};actions.append(duplicate);
  if(state.panels.filter(c=>c.species===p.id).length>1){const remove=document.createElement('button');remove.type='button';remove.textContent='Remove';remove.setAttribute('aria-label',`Remove ${meta.name} ${p.panelKey}`);remove.onclick=()=>{state.panels=state.panels.filter(c=>c!==config);load();};actions.append(remove);}label.append(actions);
  if(!p.error){const coords=document.createElement('div');coords.className='coordinate';coords.textContent=`${p.chrom}:${fmt(Math.max(1,p.from))}–${fmt(p.to)} (${p.strand})`;label.append(coords);if(p.tssTranscript){const anchor=document.createElement('div');anchor.className='coordinate';anchor.textContent=`TSS: ${p.tssTranscript} · ${p.tssMethod.replaceAll('_',' ')}`;label.append(anchor);}}
  row.append(label);const plot=document.createElement('div');plot.className='track-plot';row.append(plot);$('tracks').append(row);
  if(p.error){plot.classList.add('track-error');plot.textContent=p.error;continue;}
  const modelRows=$('models').checked?[...(p.models||[])].sort((a,b)=>(b.gname.toUpperCase()===d.gene)-(a.gname.toUpperCase()===d.gene)):[];
  const chipRows=(p.chipTracks||[]).filter(t=>t.available&&$('chip'+t.name).checked);
  const ccre=$('ccre').checked?p.ccre:null,ccreTop=94+chipRows.length*28;
  const ccreRegions=ccre?ccre.regions.filter(r=>$('ccreType').value==='all'||r.category===$('ccreType').value):[];
  const te=$('teRepeats').checked?p.teRepeats:null,teTop=ccreTop+(ccre?32:0);
  const teType=$('teType').value,teMatch=c=>teType==='all'||(teType==='other'?['Retroposon','RC'].includes(c):c===teType);
  const teClasses=te?te.classes.filter(c=>teMatch(c.name)):[],teRegions=te?te.regions.filter(r=>teMatch(r.repeatClass)):[];
  const loopTrack=$('microc').checked?p.microc:null,loopTop=teTop+(te?20+teClasses.length*7:0);
  const arcTracks=[loopTrack,$('abc').checked?p.abc:null].filter(Boolean);
  const hic=$('humanHic').checked?p.hic:null,hicTop=loopTop+arcTracks.length*140,hicDepth=700;
  const geneTop=hicTop+(hic?hicDepth+75:0);
  const height=Math.max(116,geneTop+modelRows.length*17);
  const svg=svgEl('svg',{viewBox:`0 0 ${width} ${height}`,preserveAspectRatio:'none',role:'img','aria-label':`${meta.name} ${d.gene} ATAC signal, gene annotations, ATAC peaks, ChIP-seq regions, candidate cis-regulatory elements TE repeats, Micro-C loops, ABC links and human Hi-C heatmap`});svg.style.height=height+'px';plot.append(svg);
  const unit=p.signalUnit||'CPM';
  const values=p.signal.map(v=>v.value).filter(Number.isFinite),autoMax=Math.max(.001,...values),autoMin=Math.min(0,...values);
  const shared=$('scale').value==='shared';
  let max=shared?globalMax:(config.max??autoMax),min=shared?globalMin:(config.min??autoMin);
  // Keep a one-sided manual bound valid when automatic bounds change after navigation.
  if(max<=min){if(config.max===null)max=min+Math.max(.001,Math.abs(min)*.1);else min=max-Math.max(.001,Math.abs(max)*.1);}
  const range=document.createElement('form');range.className='signal-range';range.innerHTML=`<label>Min ${esc(unit)}<input type="number" step="any" aria-label="Minimum signal ${p.panelKey}" placeholder="${min.toPrecision(3)}" value="${config.min??''}"></label><label>Max ${esc(unit)}<input type="number" step="any" aria-label="Maximum signal ${p.panelKey}" placeholder="${max.toPrecision(3)}" value="${config.max??''}"></label><div><button type="submit">Set</button><button type="button" class="auto-range">Auto</button></div>`;
  const [minInput,maxInput]=range.querySelectorAll('input');
  range.onsubmit=e=>{e.preventDefault();const lo=minInput.value===''?null:Number(minInput.value),hi=maxInput.value===''?null:Number(maxInput.value);maxInput.setCustomValidity('');if((lo!==null&&!Number.isFinite(lo))||(hi!==null&&!Number.isFinite(hi))||(lo!==null&&hi!==null&&hi<=lo)){maxInput.setCustomValidity('Maximum must be greater than minimum.');maxInput.reportValidity();return;}config.min=lo;config.max=hi;$('scale').value='independent';render();};
  minInput.oninput=maxInput.oninput=()=>maxInput.setCustomValidity('');
  range.querySelector('.auto-range').onclick=()=>{config.min=null;config.max=null;render();};label.append(range);
  const y=v=>67-(Math.max(min,Math.min(max,v))-min)/(max-min)*53;

  for(let i=0;i<=4;i++)add(svg,'line',{x1:x(-d.window+i*d.window/2),x2:x(-d.window+i*d.window/2),y1:2,y2:height-4,stroke:'#f0f3f6'});
  add(svg,'line',{x1:pad,x2:width-pad,y1:y(0),y2:y(0),stroke:'#dee2e6'});
  let segment=[];const draw=()=>{if(!segment.length)return;const pts=segment.map(v=>`${x(v.rel).toFixed(2)},${y(v.value).toFixed(2)}`).join(' ');add(svg,'polygon',{points:`${x(segment[0].rel)},${y(0)} ${pts} ${x(segment.at(-1).rel)},${y(0)}`,fill:meta.color,'fill-opacity':.8});segment=[];};
  for(const v of p.signal){if(v.value===null)draw();else segment.push(v);}draw();
  add(svg,'text',{x:width-15,y:10,'text-anchor':'end',fill:'#6c757d','font-size':8},`${min.toPrecision(3)}–${max.toPrecision(3)} ${unit}${shared?' · shared':''}`);
  if(!values.length)add(svg,'text',{x:width/2,y:40,'text-anchor':'middle',fill:'#958c7f','font-size':12},'Outside this assembly sequence');
  const tssX=x(-d.pan);if(tssX>=pad&&tssX<=width-pad)add(svg,'line',{x1:tssX,x2:tssX,y1:0,y2:height-2,stroke:'#005aa0','stroke-dasharray':'3 3','stroke-width':1});
  const sgn=p.strand==='-'?-1:1,toRel=pos=>(pos-p.centre)*sgn;
  function rect(a,b,yy,h,color){let a1=x(toRel(a)),b1=x(toRel(b)),left=Math.max(pad,Math.min(a1,b1)),right=Math.min(width-pad,Math.max(a1,b1));if(right>=left)return add(svg,'rect',{x:left,y:yy,width:Math.max(1,right-left),height:h,fill:color,rx:.6});}
  if($('peaks').checked){for(const peak of p.peaks||[])rect(peak.start,peak.end,76,4,meta.color);add(svg,'text',{x:width-15,y:87,'text-anchor':'end',fill:'#6c757d','font-size':7},p.peakAvailable?`${p.peaks.length} peaks`:(p.peakWarning||'Peak file unavailable'));}
  chipRows.forEach((track,i)=>{const yy=94+i*28;
   add(svg,'text',{x:pad,y:yy,fill:track.color,'font-size':10},`${track.name} · ${track.available?track.regions.length+' regions':'file unavailable'}`);
   for(const region of track.regions){const bar=rect(region.start,region.end,yy+5,6,track.color);if(bar)bar.chipInfo=`${meta.name} ${track.name} · ${region.name}\n${p.chrom}:${fmt(region.start)}–${fmt(region.end)} (1-based)\n${track.assembly} · replicated peaks`;}
  });
  if(ccre){
   add(svg,'text',{x:pad,y:ccreTop,fill:'#42566a','font-size':10},`ENCODE4 cCREs · all tissues · ${ccre.available?ccreRegions.length+' regions':'file unavailable'}`);
   for(const region of ccreRegions){const bar=rect(region.start,region.end,ccreTop+6,8,region.color);if(bar)bar.chipInfo=`${region.name} · ${region.category}\n${p.chrom}:${fmt(region.start)}–${fmt(region.end)} (1-based)\nENCODE4 · hg38 · all tissues`;}
  }
  if(te){
   add(svg,'text',{x:pad,y:teTop,fill:'#42566a','font-size':10},`TE repeats${teType==='all'?'':' · '+$('teType').selectedOptions[0].textContent} · ${te.available?teRegions.length+' regions':'index unavailable'}`);
   const classMap=new Map(teClasses.map((c,i)=>[c.name,{...c,row:i}]));
   for(const region of teRegions){const c=classMap.get(region.repeatClass);if(!c)continue;const bar=rect(region.start,region.end,teTop+6+c.row*7,5,c.color);if(bar)bar.chipInfo=`NMR TE · ${region.name}\n${p.chrom}:${fmt(region.start)}–${fmt(region.end)} (1-based)\nStrand: ${region.strand} · BED score: ${region.score}\n${te.assembly}`;}
  }
  arcTracks.forEach((loopTrack,trackIndex)=>{
   const top=loopTop+trackIndex*140,isABC=loopTrack.kind==='abc';
   const loops=loopTrack.regions||[],base=top+104,color=isABC?'#9b4c91':'#1864ab';
   add(svg,'text',{x:pad,y:top+4,fill:color,'font-size':10},`${loopTrack.label||'Micro-C loops'} · ${isABC?'enhancer → gene':p.sample} · ${loops.length} of ${loopTrack.total} ${isABC?'links':'loops'}${loopTrack.total>loops.length?' (highest scores)':''}`);
   add(svg,'line',{x1:pad,x2:width-pad,y1:base,y2:base,stroke:'#cdddeb'});
   if(!loops.length)add(svg,'text',{x:pad,y:top+42,fill:'#6c757d','font-size':10},loopTrack.available?'No loop anchors in this window':'Loop file unavailable');
   for(const loop of loops){
    const mid1=(loop.start1+loop.end1)/2,mid2=(loop.start2+loop.end2)/2;
    const raw1=x(toRel(mid1)),raw2=x(toRel(mid2));
    const local1=loop.chrom1===p.chrom,local2=loop.chrom2===p.chrom;
    const inside1=local1&&mid1>=p.from&&mid1<=p.to,inside2=local2&&mid2>=p.from&&mid2<=p.to;
    const clamp=v=>Math.max(pad,Math.min(width-pad,v));
    const a=local1?clamp(raw1):(inside2?clamp(raw2):pad),b=local2?clamp(raw2):(inside1?clamp(raw1):width-pad);
    const arcHeight=Math.max(15,Math.min(76,Math.abs(b-a)*.22));
    const shape=add(svg,'path',{d:`M ${a} ${base} Q ${(a+b)/2} ${base-2*arcHeight} ${b} ${base}`,fill:'none',stroke:color,'stroke-width':1.2,'stroke-opacity':.55,...(!inside1||!inside2?{'stroke-dasharray':'4 3'}:{})});
    shape.chipInfo=`${p.sample} · ${loop.name}\n${loop.chrom1}:${fmt(loop.start1)}–${fmt(loop.end1)} ↔ ${loop.chrom2}:${fmt(loop.start2)}–${fmt(loop.end2)} (1-based)\nlog2 O/E: ${loop.log2OE} · log2 enrichment: ${loop.log2Enrichment}\nCalled in ${loop.windowsCalled}/${loop.windowsCovering} covering windows${!inside1||!inside2?'\nDashed: an anchor is outside this view':''}`;
    if(isABC)shape.chipInfo=`Hypothalamus ABC prediction · TEST\nTarget gene: ${loop.targetGene} · ABC score: ${loop.score}\nEnhancer: ${loop.chrom1}:${fmt(loop.start1)}–${fmt(loop.end1)} (1-based)\nGene anchor: ${loop.chrom2}:${fmt(loop.start2)} (1-based, source startCodon point)\nSource target coordinates: ${loop.sourceTargetStart}–${loop.sourceTargetEnd}\n${loop.name}${!inside1||!inside2?'\nDashed: an anchor is outside this view':''}`;
    if(local1){const bar=rect(loop.start1,loop.end1,base-2,4,color);if(bar)bar.chipInfo=shape.chipInfo;}
    if(local2){const bar=rect(loop.start2,loop.end2,base-2,4,color);if(bar)bar.chipInfo=shape.chipInfo;}
   }
   add(svg,'text',{x:pad,y:top+117,fill:'#6c757d','font-size':8},'Arcs connect anchors · dashed: off-screen anchor · height reflects distance, not contact strength');
  });
  if(hic){
   add(svg,'text',{x:pad,y:hicTop+5,fill:'#9b2828','font-size':10},`${hic.label} · hg38 · ${hic.available?sizeLabel(hic.resolution)+' bins · '+hic.normalization:hic.error}`);
   if(hic.available){
    const top=hicTop+25,depth=hicDepth,span=2*d.window;
    const point=(a,b)=>`${((x(toRel(a))+x(toRel(b)))/2).toFixed(2)},${(top+Math.abs(b-a)/span*depth).toFixed(2)}`;
    for(const [i,j,value] of hic.cells){
     const a=Math.max(p.from,hic.origin+i*hic.resolution+1),b=Math.min(p.to,hic.chromSize,hic.origin+(i+1)*hic.resolution);
     const c=Math.max(p.from,hic.origin+j*hic.resolution+1),e=Math.min(p.to,hic.chromSize,hic.origin+(j+1)*hic.resolution);
     if(b<a||e<c)continue;
     const t=value===null?0:Math.min(1,Math.log1p(Math.max(0,value))/Math.log1p(hic.colorMax));
     const color=value===null?'#d9dfe5':`rgb(255,${Math.round(255-210*t)},${Math.round(255-215*t)})`;
     const cell=add(svg,'polygon',{points:i===j?[point(a,a),point(b,b),point(a,b)].join(' '):[point(a,c),point(b,c),point(b,e),point(a,e)].join(' '),fill:color});
     cell.chipInfo=`HepG2 Hi-C · ${hic.accession} · hg38\n${p.chrom}:${fmt(a)}–${fmt(b)} ↔ ${p.chrom}:${fmt(c)}–${fmt(e)}\n${hic.normalization}: ${value===null?'Masked / unavailable':value.toPrecision(4)} · ${sizeLabel(hic.resolution)} bins`;
    }
    add(svg,'text',{x:pad,y:hicTop+hicDepth+45,fill:'#6c757d','font-size':9},`White → red: 0 → ${hic.colorMax.toPrecision(3)} ${hic.normalization==='NONE'?'raw contacts':'KR-balanced contacts'} · log color scale · capped at 98th percentile · gray: masked`);
    add(svg,'text',{x:pad,y:hicTop+hicDepth+59,fill:'#6c757d','font-size':8},hic.warning||'HepG2 cell line, merged replicates · both anchors in view · depth represents genomic separation');
   }
  }
  if($('models').checked){if(p.modelCount>5)add(svg,'text',{x:width-15,y:height-1,'text-anchor':'end',fill:'#6c757d','font-size':7},`5 of ${p.modelCount} genes shown`);modelRows.forEach((m,i)=>{const yy=geneTop+i*17;rect(m.start,m.end,yy,1,'#8b99a6');const starts=m.ex_s.split(',').filter(Boolean).map(Number),ends=m.ex_e.split(',').filter(Boolean).map(Number);starts.forEach((a,j)=>rect(a+1,ends[j],yy-2,5,'#596f83'));add(svg,'text',{x:Math.min(width-45,Math.max(pad,x(toRel(m.start)))),y:yy+12,fill:'#596f83','font-size':8},`${m.gname} ${m.strand===p.strand?'→':'←'}`);});}
  attachTrackPan(svg,d,width,pad);
  svg.addEventListener('mousemove',event=>{if(document.body.classList.contains('track-panning'))return;if(event.target.chipInfo){tip.textContent=event.target.chipInfo;tip.hidden=false;tip.style.left=`${Math.max(4,Math.min(event.clientX+12,innerWidth-330))}px`;tip.style.top=`${Math.min(event.clientY+12,innerHeight-85)}px`;return;}const r=svg.getBoundingClientRect(),frac=Math.max(0,Math.min(1,(event.clientX-r.left)/r.width));const rel=(frac*width-pad)/(width-2*pad)*2*d.window-d.window;const point=p.signal.reduce((a,b)=>Math.abs(a.rel-rel)<Math.abs(b.rel-rel)?a:b);tip.textContent=`${meta.name} · ${p.sample}\n${p.chrom}:${fmt(p.centre+sgn*point.rel)}\n${point.value===null?'Outside sequence':point.value.toFixed(3)+' '+unit} · ${sizeLabel(point.rel+d.pan)} from TSS`;tip.hidden=false;tip.style.left=`${Math.max(4,Math.min(event.clientX+12,innerWidth-310))}px`;tip.style.top=`${Math.min(event.clientY+12,innerHeight-85)}px`;});svg.addEventListener('mouseleave',()=>tip.hidden=true);
 }
}
$('searchForm').onsubmit=e=>{e.preventDefault();state.gene=$('gene').value.trim().toUpperCase();state.pan=0;load();};
let searchTimer;$('gene').oninput=()=>{clearTimeout(searchTimer);const q=$('gene').value.trim();$('genes').replaceChildren();$('annotationCount').textContent=q?'Checking annotations…':'';searchTimer=setTimeout(async()=>{try{const result=await json('/api/genes?q='+encodeURIComponent(q));if(q!==$('gene').value.trim())return;$('genes').replaceChildren(...result.choices.map(a=>{const option=new Option(a.gene,a.gene);option.label=`${a.gene} · ${a.samples} samples / ${a.species} species`;return option;}));$('annotationCount').textContent=q?annotationLabel(result.match):'';}catch{$('annotationCount').textContent='Annotation counts unavailable';}},180);};
function attachTrackPan(svg,data,width,pad){
 svg.addEventListener('pointerdown',event=>{
  if(event.button!==0||!event.isPrimary||state.loading)return;
  const startX=event.clientX,plotPixels=svg.getBoundingClientRect().width*(width-2*pad)/width;
  if(plotPixels<=0)return;
  const views=[...document.querySelectorAll('#tracks .track-plot svg,#ruler svg')].map(el=>({el,box:el.getAttribute('viewBox')}));
  let delta=0,moved=false;
  const oldStatus=$('status').textContent;
  svg.setPointerCapture(event.pointerId);tip.hidden=true;
  document.body.classList.add('track-panning');$('exportPdf').disabled=true;
  const move=e=>{
   if(e.pointerId!==event.pointerId)return;
   const dx=e.clientX-startX;
   if(!moved&&Math.abs(dx)<4)return;
   moved=true;delta=-Math.round(dx/plotPixels*data.window*2);
   for(const v of views){const box=v.box.split(' ').map(Number);box[0]=delta/(data.window*2)*(width-2*pad);v.el.setAttribute('viewBox',box.join(' '));}
   $('status').textContent=`Pan ${sizeLabel(delta)} · release to load · Esc to cancel`;
  };
  const finish=(commit)=>{
   svg.removeEventListener('pointermove',move);svg.removeEventListener('pointerup',up);svg.removeEventListener('pointercancel',cancel);svg.removeEventListener('lostpointercapture',cancel);document.removeEventListener('keydown',key);
   if(svg.hasPointerCapture(event.pointerId))svg.releasePointerCapture(event.pointerId);
   for(const v of views)v.el.setAttribute('viewBox',v.box);
   document.body.classList.remove('track-panning');$('status').textContent=oldStatus;$('exportPdf').disabled=state.loading||!state.data;
   if(commit&&moved&&delta&&state.data===data&&!state.loading){state.pan=data.pan+delta;load();}
  };
  const up=e=>{if(e.pointerId===event.pointerId){move(e);finish(true);}};
  const cancel=()=>finish(false),key=e=>{if(e.key==='Escape'){e.preventDefault();finish(false);}};
  svg.addEventListener('pointermove',move);svg.addEventListener('pointerup',up);svg.addEventListener('pointercancel',cancel);svg.addEventListener('lostpointercapture',cancel);document.addEventListener('keydown',key);
 });
}
function windowChange(span){state.window=Math.max(1,Math.round(span))/2;$('window').value=state.window*2;load();}
$('viewForm').onsubmit=e=>{e.preventDefault();state.bins=Number($('bins').value);windowChange(Number($('window').value));};
$('zoomIn').onclick=()=>windowChange(state.window);$('zoomOut').onclick=()=>windowChange(state.window*4);
$('left').onclick=()=>{state.pan-=Math.round(state.window*.5);load();};$('right').onclick=()=>{state.pan+=Math.round(state.window*.5);load();};$('reset').onclick=()=>{state.pan=0;windowChange(10000);};
for(const id of ['models','peaks','scale','chipH3K4me3','chipH3K27ac','teRepeats','teType','microc','abc','humanHic','ccre','ccreType'])$(id).onchange=render;
$('exportPdf').onclick=async()=>{
 if(!state.data||state.loading)return;
 const button=$('exportPdf');button.disabled=true;button.textContent='Creating PDF…';
 const d=state.data,serializer=new XMLSerializer();
 const payload={gene:d.gene,summary:$('summary').textContent+' · '+$('scale').selectedOptions[0].textContent,
  ruler:serializer.serializeToString($('ruler').querySelector('svg')),
  panels:d.panels.map(p=>{const row=document.querySelector(`[data-panel-key="${p.panelKey}"]`),svg=row.querySelector('svg');return {name:state.catalog.find(s=>s.id===p.id).name,sample:p.sample||'',coordinates:row.querySelector('.coordinate')?.textContent||'',error:p.error||'',svg:svg?serializer.serializeToString(svg):null};})};
 try{const response=await fetch('/api/export-pdf',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});if(!response.ok){const error=await response.json();throw Error(error.error);}
 const result=await response.json(),a=document.createElement('a');a.href=result.url;a.download=`${d.gene}_NSC.pdf`;document.body.append(a);a.click();a.remove();
 }catch(e){$('status').textContent=e.message;}finally{button.textContent='↓ Export PDF';button.disabled=state.loading||!state.data;}
};

function saveSpeciesOrder(){
 orderPanels(state.panels);try{localStorage.setItem('nscSpeciesOrder',JSON.stringify(state.catalog.map(x=>x.id)));}catch{}
 renderSpecies();render();
}
function attachSpeciesDrag(handle,row,species){
 handle.onpointerdown=e=>{
  if(e.button!==0)return;e.preventDefault();handle.setPointerCapture(e.pointerId);
  let target=null,after=false,moved=false;
  const clear=()=>document.querySelectorAll('.species-order-row').forEach(r=>r.classList.remove('drop-before','drop-after'));
  const move=event=>{
   if(Math.abs(event.clientY-e.clientY)<4&&!moved)return;moved=true;row.classList.add('dragging');clear();
   const sidebar=document.querySelector('.sidebar'),box=sidebar.getBoundingClientRect();
   if(event.clientY<box.top+40)sidebar.scrollTop-=18;
   if(event.clientY>box.bottom-40)sidebar.scrollTop+=18;
   const hit=document.elementFromPoint(event.clientX,event.clientY)?.closest('.species-order-row');
   target=hit&&hit!==row?hit:null;
   if(target){const r=target.getBoundingClientRect();after=event.clientY>r.top+r.height/2;target.classList.add(after?'drop-after':'drop-before');}
  };
  const finish=event=>{
   if(event.type==='pointerup')move(event);const id=target?.dataset.species;clear();row.classList.remove('dragging');document.removeEventListener('pointermove',move);document.removeEventListener('pointerup',finish);document.removeEventListener('pointercancel',finish);document.removeEventListener('keydown',cancel);
   if(handle.hasPointerCapture(e.pointerId))handle.releasePointerCapture(e.pointerId);
   if(event.type==='pointerup'&&moved&&id){const from=state.catalog.findIndex(s=>s.id===species.id),[item]=state.catalog.splice(from,1);const to=state.catalog.findIndex(s=>s.id===id)+(after?1:0);state.catalog.splice(to,0,item);saveSpeciesOrder();}
  };
  const cancel=event=>{if(event.key==='Escape')finish(event);};
  document.addEventListener('keydown',cancel);document.addEventListener('pointermove',move);document.addEventListener('pointerup',finish);document.addEventListener('pointercancel',finish);
 };
}
function renderSpecies(){
 $('species').replaceChildren();
 state.catalog.forEach((s,i)=>{
  const row=document.createElement('div');row.className='species-order-row';row.dataset.species=s.id;
  const grip=document.createElement('span');grip.className='species-drag-handle';grip.textContent='⠿';grip.title=`Drag to reorder ${s.name}`;grip.setAttribute('aria-hidden','true');row.append(grip);attachSpeciesDrag(grip,row,s);
  const label=document.createElement('label');label.className='species-choice';
  const check=document.createElement('input');check.type='checkbox';check.checked=state.species.includes(s.id);
  check.onchange=()=>{state.species=state.catalog.filter(x=>x.id===s.id?check.checked:state.species.includes(x.id)).map(x=>x.id);load();};
  const text=document.createElement('span');text.innerHTML=`${esc(s.name)}<small>${esc(s.latin)}</small>`;label.append(check,text);row.append(label);
  $('species').append(row);
 });
}
(async()=>{try{
 state.catalog=await json('/api/catalog');
 try{const saved=JSON.parse(localStorage.getItem('nscSpeciesOrder')||'[]');if(Array.isArray(saved))state.catalog.sort((a,b)=>(saved.includes(a.id)?saved.indexOf(a.id):saved.length)-(saved.includes(b.id)?saved.indexOf(b.id):saved.length));}catch{}
 state.species=state.catalog.map(s=>s.id);for(const s of state.catalog)state.panels.push(newPanel(s));renderSpecies();await load();
}catch(e){$('status').textContent=`Could not connect to the Python server: ${e.message}`;}})();
