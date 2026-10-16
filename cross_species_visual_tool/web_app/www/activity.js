/* Precomputed gene-region activity; one independent observation per species. */
function createActivityPanel(prefix, handle) {
 const el=id=>document.getElementById(id.replace(/^activity/,prefix)), NS='http://www.w3.org/2000/svg';
 let data=null, requested='', serial=0;
 const shortNames={Hydrochoerus_hydrochaeris:'Capybara',Eonycteris_spelaea:'Nectar bat',Myocastor_coypus:'Nutria',Balaena_mysticetus:'Bowhead',Artibeus_jamaicensis:'Fruit bat',Castor_canadensis:'Beaver',Cricetulus_griseus:'Chinese hamster',Chinchilla_lanigera:'Chinchilla',Octodon_degus:'Degu',Bramus_lutescens:'B. lutescens',Eptesicus_fuscus:'Big brown bat',Sciurus_niger:'Fox squirrel',Meriones_unguiculatus:'Gerbil',Mesocricetus_auratus:'Golden hamster',Molossus_molossus:'Mastiff bat',Erethizon_dorsatum:'Porcupine',Tamiasciurus_hudsonicus:'Red squirrel',Marmota_monax:'Woodchuck',Peromyscus_maniculatus:'Deer mouse',Microtus_pennsylvanicus:'Meadow vole',Cavia_porcellus:'Guinea pig',Sciurus_carolinensis:'Gray squirrel'};
 const fmt=x=>Number.isFinite(x)?Number(x.toPrecision(4)).toLocaleString('en-US'):'—';
 function element(tag,attrs={},text){const n=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n;}
 function row(root,values){const tr=document.createElement('tr');for(const value of values){const td=document.createElement('td');td.textContent=value;tr.append(td);}root.append(tr);}
 function render(){
  if(!data)return;
  const catalog=state.catalog.filter(s=>el('activityScope').value==='all'||state.species.includes(s.id));
  const peakSignal=data.strategy.startsWith('peak_signal'), regulatory=data.strategy.startsWith('peaks');
  const peakWindow=Number(data.strategy.match(/(?:peaks|peak_signal)(\d+)kb$/)?.[1]);
  if(prefix==='strength'){
   const scale=el('activityScale');
   scale.options[0].textContent=peakSignal?'Median CPM signal in peaks':'Median peak count';
   scale.options[1].textContent=peakSignal?'log₁₀(Median CPM signal in peaks)':'log₁₀(Median peak count)';
  }
  const axis=el('activityAxis').value, logActivity=el('activityScale').value==='log10', points=[];
  el('activityValueHeading').textContent=regulatory?(logActivity?'log₁₀(Median peak count)':'Median peak count'):peakSignal?(logActivity?'log₁₀(Median signal in peaks)':'Median signal in peaks'):logActivity?'log₁₀(Median activity)':'Median activity';
  const root=el('activityChart');root.replaceChildren();el('activityRows').replaceChildren();el('activitySamples').replaceChildren();
  el('activityGene').textContent=data.gene;
  el('activityRegion').textContent=peakSignal?`CPM signal in peaks within TSS ±${peakWindow} kb`:regulatory?`Regulatory activity: peak count within TSS ±${peakWindow} kb`:data.strategy==='body'?'Gene-body + 2 kb upstream signal':`TSS ±${{tss2kb:2,tss5kb:5,tss100kb:100,tss200kb:200,tss500kb:500}[data.strategy]} kb signal`;
  for(const sp of catalog){
   const records=data.records.filter(r=>r.species===sp.id), life=data.lifespans[sp.latin];
   if(!records.length){row(el('activityRows'),[sp.name,'—',life?.display_years||fmt(life?.years),fmt(life?.body_mass_g),'No gene match']);continue;}
   const selected=[...records].sort((a,b)=>{
    const av=Number.isFinite(a.median),bv=Number.isFinite(b.median);
    return av!==bv?(av?-1:1):av&&a.median!==b.median?b.median-a.median:a.locus.localeCompare(b.locus);
   })[0];
   for(const r of records)for(const [sample,value] of Object.entries(r.samples))row(el('activitySamples'),[sp.name,r.locus+(records.length>1?(r===selected?' (selected)':' (not selected)'):''),sample,fmt(value),r.unit,r.window_details?.[sample]?.quantified_bp??'—',r.window_details?.[sample]?.status??'—']);
   for(const r of [selected]){
    const y=logActivity?(r.median>0?Math.log10(r.median):null):r.median;
    const span=r.locus.match(/:(\d+)-(\d+):[+-]$/);
    const lengthKb=span?(Number(span[2])-Number(span[1]))/1000:NaN;
    const reason=!Number.isFinite(r.median)?(regulatory||peakSignal?'Signal or peaks unavailable / incompatible':'No measurement'):logActivity&&r.median<=0?'Nonpositive activity (log undefined)':sp.id==='Human'?(regulatory?'Excluded from fit':'Different signal scale'):sp.id==='Macaque'?'Excluded from fit':sp.id==='Balaena_mysticetus'&&el('activityExcludeBowhead').checked?'Excluded by option':!Number.isFinite(life?.[axis])?(life?.lower_bound_years?'MLS lower bound':'MLS unavailable'):data.strategy==='body'&&(!Number.isFinite(lengthKb)||lengthKb<=0)?'Gene length unavailable':'Included';
    row(el('activityRows'),[sp.name,fmt(y)+(logActivity?'':regulatory?' peaks':r.unit==='CPM_signal_bp'?' CPM·bp':' scaled·bp'),life?.display_years||fmt(life?.years),fmt(life?.body_mass_g),reason+(records.length>1?' · highest-activity locus':'')]);
    if(reason==='Included')points.push({x:life[axis],y,rawActivity:r.median,sp,label:shortNames[sp.id]||sp.id,unit:r.unit,lengthKb,logMass:life.log10_body_mass_g});
   }
  }
  el('activityStatus').textContent=`${points.length} species included in the correlation. For multiple loci, the highest species-median activity is used. Values use the selected activity region and do not change when the browser window is panned or zoomed.`;
  if(el('activityExcludeBowhead').checked)el('activityStatus').textContent+=' Bowhead whale is excluded from the plot and correlation; its activity remains in the table.';
  if(regulatory)el('activityStatus').textContent+=' Regulatory activity is the median number of overlapping peaks; peak calling and sequencing depth can affect this score. No covariates are included.';
  if(peakSignal)el('activityStatus').textContent+=' Signal is integrated only over the union of sample peaks inside the TSS window (CPM·bp), without distance weighting or covariates. Peak calling and sequencing depth can affect this score.';
  if(logActivity)el('activityStatus').textContent+=' Plot and correlation use log₁₀ of the species median; nonpositive activity is excluded without a pseudocount.';
  if(!regulatory&&points.some(p=>p.unit!=='CPM_signal_bp'))el('activityStatus').textContent+=' Mixed source scales: human uses scaled·bp; other species use CPM·bp. No cross-species rescaling is applied.';
  if(!points.length){root.textContent='No comparable gene activity values with MLS for this selection.';return;}
  const svg=element('svg',{viewBox:'0 0 620 390',role:'img','aria-label':`${data.gene} ${logActivity?'log10 median':'median'} ${peakSignal?'signal in peaks':regulatory?'regulatory activity':'gene activity'} versus ${axis==='years'?'MLS':'log10 MLS'}`});
  const left=90,right=595,top=35,bottom=315;
  const extent=values=>{const min=Math.min(...values),max=Math.max(...values),pad=(max-min||Math.abs(max)||1)*.15;return [min-pad,max+pad];};
  let [xmin,xmax]=extent(points.map(p=>p.x)),[ymin,ymax]=extent(points.map(p=>p.y));
  if(axis==='years')xmin=Math.max(0,xmin);if(!logActivity)ymin=Math.max(0,ymin);
  const X=v=>left+(v-xmin)/(xmax-xmin)*(right-left),Y=v=>bottom-(v-ymin)/(ymax-ymin)*(bottom-top);
  for(let i=0;i<=4;i++){
   const x=xmin+(xmax-xmin)*i/4,y=ymin+(ymax-ymin)*i/4;
   svg.append(element('line',{x1:left,x2:right,y1:Y(y),y2:Y(y),stroke:'#e2e8ee'}),element('text',{x:left-9,y:Y(y)+4,'text-anchor':'end','font-size':11},fmt(y)),element('text',{x:X(x),y:bottom+23,'text-anchor':'middle','font-size':11},fmt(x)));
  }
  svg.append(element('text',{x:(left+right)/2,y:370,'text-anchor':'middle','font-size':13},axis==='years'?'Maximum lifespan (years)':'log₁₀(Maximum lifespan in years)'),element('text',{transform:'translate(18 175) rotate(-90)','text-anchor':'middle','font-size':13},regulatory?(logActivity?'log₁₀(Median peak count)':'Regulatory activity (median peak count)'):peakSignal?(logActivity?'log₁₀(Median signal in peaks)':'Median signal in peaks (CPM·bp)'):logActivity?'log₁₀(Median gene activity)':points.some(p=>p.unit!=='CPM_signal_bp')?'Median gene activity (source units·bp)':'Median gene activity (CPM·bp)'));
  const adjusted=data.strategy==='body';
  const fit=adjusted?NSCStatistics.logLengthAdjustedFit(points):NSCStatistics.linearFit(points);
  if(adjusted)el('activityStatus').textContent+=' Gene-body fit adjusts only for log₁₀(annotated gene length in kb), without body mass. Points show activity on the selected scale; the line holds log₁₀(gene length) at its species mean.';
  if(fit){
   let a=Math.min(...points.map(p=>p.x)),b=Math.max(...points.map(p=>p.x));
   if(fit.slope!==0){const ends=[(ymin-fit.intercept)/fit.slope,(ymax-fit.intercept)/fit.slope].sort((a,b)=>a-b);a=Math.max(a,ends[0]);b=Math.min(b,ends[1]);}
   if(a<=b)svg.append(element('line',{x1:X(a),x2:X(b),y1:Y(fit.slope*a+fit.intercept),y2:Y(fit.slope*b+fit.intercept),stroke:'#005aa0','stroke-width':2}));
  }
  const labels=[];
  for(const p of points){
   const cx=X(p.x),cy=Y(p.y),width=p.label.length*7+8;
   let ly=cy-10;
   while(labels.some(q=>Math.abs(q.x-cx)<(q.width+width)/2&&Math.abs(q.y-ly)<14))ly-=15;
   labels.push({x:cx,y:ly,width});
   if(ly<cy-10)svg.append(element('line',{x1:cx,y1:cy-6,x2:cx,y2:ly+3,stroke:p.sp.color,'stroke-width':.7}));
   const c=element('circle',{cx,cy,r:5,fill:p.sp.color,tabindex:0});c.append(element('title',{},`${p.sp.name}: ${logActivity?'log₁₀(median activity) '+fmt(p.y)+'; ':''}median activity ${fmt(p.rawActivity)} ${regulatory?'peaks':p.unit==='CPM_signal_bp'?'CPM·bp':'scaled·bp'}; ${axis} ${fmt(p.x)}`));svg.append(c,element('text',{x:cx,y:ly,'text-anchor':'middle','font-size':11},p.label));
  }
  const stats=document.createElement('p');stats.className='omics-stats';
  stats.textContent=fit?`Pearson r = ${fit.r.toFixed(3)} · p = ${fit.p===null?'unavailable':fit.p===0?'< 1e-300':fit.p<.001?fit.p.toExponential(2):fit.p.toFixed(4)} · R² = ${(fit.r**2).toFixed(3)} · n = ${points.length} species · y = ${fmt(fit.slope)}x ${fit.intercept<0?'−':'+'} ${fmt(Math.abs(fit.intercept))}`:'A fit requires at least 3 species with variable lifespan and activity.';
  if(adjusted)stats.textContent=fit?`Log-length-adjusted partial r = ${fit.r.toFixed(3)} · p = ${fit.p===null?'unavailable':fit.p<.001?fit.p.toExponential(2):fit.p.toFixed(4)} · partial R² = ${(fit.r**2).toFixed(3)} · n = ${points.length} species · df = ${fit.df} · MLS coefficient = ${fmt(fit.slope)} · length coefficient = ${fmt(fit.lengthSlope)} per log₁₀(kb)`:'Adjusted fit requires at least 5 species with complete MLS, positive gene length, variable covariates and non-collinear predictors.';
  stats.title=adjusted?'Model: activity = intercept + MLS-axis coefficient × MLS-axis value + length coefficient × log10(annotated gene length in kb). Body mass is not included. Partial correlation controls for log10 gene length; two-sided t test, df = n − 3. No phylogenetic or multiple-testing adjustment.':'Two-sided Pearson test, species-level medians, unadjusted for phylogeny and multiple testing.';
  stats.title+=` Response: ${logActivity?'log10(species median activity)':'species median activity'}.`;
  root.append(svg,stats);
 }
 async function show(gene){
  gene=gene.trim().toUpperCase();let strategy=el('activityStrategy').value;
  if(prefix==='strength'&&el('activityMeasure').value==='signal')strategy=strategy.replace(/^peaks/,'peak_signal');
  const key=gene+':'+strategy;if(key===requested&&data){render();return;}if(key===requested)return;
  requested=key;const token=++serial;data=null;el('activityStatus').textContent='Loading gene activity…';
  for(const id of ['activityChart','activityRows','activitySamples'])el(id).replaceChildren();
  try{const response=await fetch('/api/gene-activity?gene='+encodeURIComponent(gene)+'&strategy='+strategy);const result=await response.json();if(token!==serial)return;if(!response.ok)throw Error(result.error||'Unable to load gene activity');data=result;render();}
  catch(error){if(token===serial){requested='';el('activityStatus').textContent=error.message;}}
 }
 el('activityStrategy').onchange=()=>show(state.gene);
 if(prefix==='strength')el('activityMeasure').onchange=()=>show(state.gene);
 el('activityAxis').onchange=el('activityScope').onchange=el('activityScale').onchange=el('activityExcludeBowhead').onchange=render;
 window[handle]={show};show(state.gene);
}
createActivityPanel('activity','NSCActivity');
createActivityPanel('strength','NSCStrength');
