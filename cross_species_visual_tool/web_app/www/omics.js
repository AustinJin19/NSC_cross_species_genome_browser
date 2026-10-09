(() => {
 const el=id=>document.getElementById(id), NS='http://www.w3.org/2000/svg';
 let data=null, requested='', serial=0;
 const format=x=>x===null||x===undefined?'Missing':Number(x).toLocaleString('en-US',{maximumSignificantDigits:5});
 const tick=v=>Math.abs(v)>=10000||(v!==0&&Math.abs(v)<.001)?v.toExponential(2):Number(v.toPrecision(3)).toLocaleString();
 const names={'Heterocephalus glaber':'Naked mole-rat (NMR)','Nannospalax galili':'Blind mole-rat (BMR)','Fukomys damarensis':'Damaraland mole-rat (DMR)','Homo sapiens':'Human','Mus musculus':'Mouse','Rattus norvegicus':'Rat','Macaca mulatta':'Rhesus macaque','Sylvilagus floridanus':'Eastern cottontail','Acomys cahirinus':'African spiny mouse'};
 const shortNames={'Heterocephalus glaber':'NMR','Nannospalax galili':'BMR','Fukomys damarensis':'DMR','Homo sapiens':'Human','Mus musculus':'Mouse','Rattus norvegicus':'Rat','Macaca mulatta':'Rhesus macaque','Sylvilagus floridanus':'Eastern cottontail','Acomys cahirinus':'African spiny mouse'};
 const palette={rna:'#005aa0',protein:'#967338'};
 function node(tag,attrs={},text){const n=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n;}
 function append(svg,tag,attrs,text){const n=node(tag,attrs,text);svg.append(n);return n;}
 function median(values){const v=[...values].sort((a,b)=>a-b),i=Math.floor(v.length/2);return v.length%2?v[i]:(v[i-1]+v[i])/2;}
 function selectedSpecies(){return state.catalog.filter(s=>state.species.includes(s.id)).map(s=>s.latin);}
 function render(){
  if(!data)return;
  const selected=el('omicsScope').value==='selected';
  const species=selected?selectedSpecies():[...new Set(data.catalog.map(r=>r.species))].sort();
  const rows=data.records.filter(r=>species.includes(r.species));
  el('omicsGene').textContent=data.gene;el('omicsExport').disabled=!data.records.length;
  el('omicsStatus').textContent=data.stale.length?`Source files changed: ${data.stale.join(', ')}. Rebuild the liver index to update this view.`:!species.length?'Select at least one ATAC species.':!rows.length?`No liver records for ${data.gene} in the selected species.`:'';
  for(const assay of ['rna','protein']){
   const metric=el(assay==='rna'?'rnaMetric':'omicsMetric').value;
   const records=rows.filter(r=>r.assay===assay), valid=records.filter(r=>Number.isFinite(r[metric]));
   const normOn=[...new Set(records.map(r=>r.norm_on).filter(Boolean))].join(', ');
   const unit=metric==='raw'?(assay==='rna'?'TPM':'Protein quant (source units unspecified)'):metric==='log'?(assay==='rna'?'log.tpm (source values)':'log.quant (source values)'):`quant.norm (normalized ${normOn||(assay==='rna'?'log.tpm':'log.quant')})`;
   el(assay+'Unit').textContent=unit+' · '+[...new Set(records.map(r=>r.workbook))].join(' + ');
   el(assay+'Summary').textContent=`${new Set(valid.map(r=>r.species)).size} species · ${new Set(valid.map(r=>r.species+'|'+r.sample)).size} samples`;
   const root=el(assay+'Chart');root.replaceChildren();
   if(!species.length){root.textContent='No species selected';continue;}
   if(!valid.length){root.textContent=`No ${metric==='raw'?'raw':metric} ${assay==='rna'?'RNA':'protein'} measurements for ${data.gene} in Liver.`;continue;}
   const axis=el('omicsAxis').value,xLabel=axis==='years'?'Maximum lifespan (years)':'log₁₀(Maximum lifespan in years)';
   const points=species.map(sp=>{const group=valid.filter(r=>r.species===sp),life=data.lifespans[sp];return group.length&&life&&Number.isFinite(life[axis])?{species:sp,x:life[axis],y:median(group.map(r=>r[metric])),n:group.length,life}:null;}).filter(Boolean);
   const missing=new Set(valid.map(r=>r.species)).size-points.length;
   el(assay+'Summary').textContent=`${points.length} species plotted · ${valid.length} records`;
   if(!points.length){root.textContent='No species have both this measurement and a verified lifespan.';continue;}
   const width=560,height=360,left=65,right=540,top=22,bottom=298;
   const svg=node('svg',{viewBox:`0 0 ${width} ${height}`,role:'img','aria-label':`${data.gene} Liver ${assay}: species median abundance versus ${xLabel}`});svg.style.width='100%';root.append(svg);
   function extent(values){let lo=Math.min(...values),hi=Math.max(...values),pad=(hi-lo||Math.abs(lo)||1)*.09;return[lo-pad,hi+pad];}
   let [xmin,xmax]=extent(points.map(p=>p.x));if(axis==='years')xmin=Math.max(0,xmin);
   const [ymin,ymax]=extent(points.map(p=>p.y));
   const x=v=>left+(v-xmin)/(xmax-xmin)*(right-left),y=v=>bottom-(v-ymin)/(ymax-ymin)*(bottom-top);
   for(let i=0;i<=4;i++){
    const xv=xmin+(xmax-xmin)*i/4,yv=ymin+(ymax-ymin)*i/4;
    append(svg,'line',{x1:left,x2:right,y1:y(yv),y2:y(yv),stroke:'#edf0f3'});
    append(svg,'text',{x:left-9,y:y(yv)+4,'text-anchor':'end','font-size':11,fill:'#6c757d'},tick(yv));
    append(svg,'text',{x:x(xv),y:bottom+19,'text-anchor':'middle','font-size':11,fill:'#6c757d'},tick(xv));
   }
   append(svg,'line',{x1:left,x2:left,y1:top,y2:bottom,stroke:'#b7c4d0'});append(svg,'line',{x1:left,x2:right,y1:bottom,y2:bottom,stroke:'#b7c4d0'});
   append(svg,'text',{x:(left+right)/2,y:height-10,'text-anchor':'middle','font-size':12,fill:'#33465a'},xLabel);
   append(svg,'text',{transform:`translate(15 ${(top+bottom)/2}) rotate(-90)`,'text-anchor':'middle','font-size':11,fill:'#33465a'},'Species median abundance');
   const fit=NSCStatistics.linearFit(points);
   if(fit){
    // Clip the fitted segment to the plot range without altering its coefficients.
    let a=Math.min(...points.map(p=>p.x)),b=Math.max(...points.map(p=>p.x));
    if(fit.slope!==0){const bounds=[(ymin-fit.intercept)/fit.slope,(ymax-fit.intercept)/fit.slope].sort((a,b)=>a-b);a=Math.max(a,bounds[0]);b=Math.min(b,bounds[1]);}
    if(a<=b)append(svg,'line',{x1:x(a),x2:x(b),y1:y(fit.intercept+fit.slope*a),y2:y(fit.intercept+fit.slope*b),stroke:palette[assay],'stroke-width':2,'stroke-opacity':.7});
   }
   for(const p of points){const circle=append(svg,'circle',{cx:x(p.x),cy:y(p.y),r:5,fill:palette[assay],'fill-opacity':.75,tabindex:0});
    const tip=`${names[p.species]?names[p.species]+'\n':''}${p.species}\nMLS: ${format(p.life.years)} years; log₁₀(MLS): ${format(p.life.log10_years)}\nMedian ${unit}: ${format(p.y)}\n${p.n} liver records\nRNA/protein source: ${[...new Set(valid.filter(r=>r.species===p.species).map(r=>r.workbook))].join(', ')}\n${p.life.source}`;
    append(circle,'title',{},tip);circle.setAttribute('aria-label',tip);
    if(shortNames[p.species])append(svg,'text',{x:x(p.x),y:y(p.y)-10,'text-anchor':'middle','font-size':10,fill:palette[assay],'paint-order':'stroke',stroke:'#fff','stroke-width':3,'stroke-linejoin':'round','pointer-events':'none'},shortNames[p.species]);
    circle.onmouseenter=circle.onfocus=()=>{detail.textContent=tip;};circle.onmouseleave=circle.onblur=()=>{detail.textContent='Hover or focus a point to inspect its species and values.';};
   }
   const stats=document.createElement('p');stats.className='omics-stats';
   stats.title='Two-sided Pearson correlation test, using species medians; df = n − 2. Unadjusted p-value.';
   const pText=fit?.p===null?'unavailable':fit?.p===0?'< 1e-300':fit?.p<.001?fit.p.toExponential(2):fit?.p.toFixed(4);
   stats.textContent=fit?`Pearson r = ${fit.r.toFixed(3)} · p ${fit.p===0?'':'= '}${pText} · R² = ${(fit.r*fit.r).toFixed(3)} · n = ${points.length} species · y = ${format(fit.slope)}x ${fit.intercept<0?'−':'+'} ${format(Math.abs(fit.intercept))}`:`n = ${points.length} species · Linear correlation requires at least 3 species and variation on both axes.`;
   root.append(stats);
   if(missing){const notice=document.createElement('p');notice.textContent=`${missing} species with abundance excluded: no verified MLS (${species.filter(sp=>valid.some(r=>r.species===sp)&&!Number.isFinite(data.lifespans[sp]?.[axis])).join(', ')}).`;root.append(notice);}
   if(selected){const absent=species.filter(sp=>!valid.some(r=>r.species===sp));if(absent.length){const note=document.createElement('p');note.textContent='No measurement in this source: '+absent.map(sp=>names[sp]||sp).join(', ')+'.';root.append(note);}}
   const detail=document.createElement('p');detail.className='omics-point-detail';detail.textContent='Hover or focus a point to inspect its species and values.';root.append(detail);
  }
  const tbody=el('omicsRows');tbody.replaceChildren();
  for(const r of rows){const tr=document.createElement('tr');for(const v of [r.assay==='rna'?'RNA':'Protein',r.species,r.sample,format(r.raw),format(r.log),format(r.normalized),r.workbook,r.source_row]){const td=document.createElement('td');td.textContent=v;tr.append(td);}tbody.append(tr);}
 }
 async function show(gene){
  gene=String(gene).trim().toUpperCase();if(requested===gene){if(data)render();return;}
  requested=gene;data=null;const current=++serial;el('omicsGene').textContent=gene;el('omicsStatus').textContent='Loading liver measurements…';el('omicsExport').disabled=true;
  for(const id of ['rnaChart','proteinChart','omicsRows'])el(id).replaceChildren();
  try{const response=await fetch('/api/omics?gene='+encodeURIComponent(gene)),result=await response.json();if(current!==serial)return;if(!response.ok)throw Error(result.error);data=result;render();}
  catch(e){if(current===serial){requested='';el('omicsStatus').textContent=e.message;}}
 }
 el('rnaMetric').onchange=el('omicsScope').onchange=el('omicsMetric').onchange=el('omicsAxis').onchange=render;
 el('omicsExport').onclick=()=>{if(!data)return;const a=document.createElement('a');a.href='/api/omics-export?gene='+encodeURIComponent(data.gene);a.download=data.gene+'_liver_omics.csv';document.body.append(a);a.click();a.remove();};
 window.NSCOmics={show};show(state.gene);
})();
