import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const out=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const data=JSON.parse(await fs.readFile(path.join(out,'matrix_data.json'),'utf8'));
const wb=Workbook.create();
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const species=data.metadata.species, last=col(3+species.length), firstRow=8, lastRow=70, sourceFirst=80, sourceLast=142;
const headerColor='#183F62', ink='#203246', pale='#EDF3F8';
const configs=[
 {key:'rna',name:'RNA expression',title:'Liver RNA expression ranks',unit:'log.tpm',
  coverage:'46 of 63 labels have measurements. Liver RNA sources cover 23 of the 32 species.',
  source:'Sources: CrossSpecies_Count.xlsx (liver); rnaLongTPM.xlsx (mouse liver supplement).',
  definition:'Values: species median of source log.tpm. No additional logarithm is applied.',format:'0.000'},
 {key:'atac',name:'TSS 5kb activity',title:'TSS ±5 kb gene activity ranks',unit:'CPM·bp',
  coverage:'52 of 63 labels have eligible measurements across the 32 ATAC species.',
  source:'Sources: current browser gene-locus index and sample CPM bigWigs. Calculated 29 October 2026.',
  definition:'Values: species median integrated CPM over TSS ±5 kb (10 kb total). No distance weighting.',format:'#,##0.0'}
];

for(const cfg of configs){
 const sh=wb.worksheets.add(cfg.name), rows=data[cfg.key];
 sh.showGridLines=false;sh.tabColor=cfg.key==='rna'?'#183F62':'#375F49';
 const all=sh.getRange(`A1:${last}157`);
 all.format.font={name:'Arial',size:10,color:ink};all.format.rowHeight=19;all.format.verticalAlignment='center';
 sh.getRange('A1:A157').format.columnWidth=21;
 sh.getRange('B1:B157').format.columnWidth=15;
 sh.getRange('C1:C157').format.columnWidth=18;
 sh.getRange(`D1:${last}157`).format.columnWidth=17;
 for(const r of [2,74])sh.getRange(`A${r}:${last}${r}`).format.font={name:'Arial',size:14,bold:true,color:headerColor};
 sh.getRange('A2').values=[[cfg.title]];
 sh.getRange('A3').values=[['Rank 1 = highest among measured screenshot labels within each species. Exact ties use average ranks. NA = unavailable.']];
 sh.getRange('A4').values=[['63 labels; 32 species. Human and bowhead whale excluded; macaque included where data are available.']];
 sh.getRange('A5').values=[['Rows ordered by mean rank across measured species (lower first). Coverage varies; compare the Species measured column.']];
 sh.getRange('A6').values=[[cfg.coverage]];
 sh.getRange(`A2:${last}2`).format.borders={bottom:{style:'thin',color:'#9FB1C0'}};
 sh.getRange(`A3:${last}6`).format.font={name:'Arial',size:10,color:'#536579'};
 sh.getRange(`A7:${last}7`).values=[['Gene / label','Mean rank','Species measured',...species.map(s=>s.common)]];
 sh.getRange(`A8:A70`).values=rows.map(r=>[r.gene]);
 sh.getRange(`D80:${last}142`).values=rows.map(r=>r.values.map(v=>v===null?'NA':v));
 sh.getRange('A80:A142').values=rows.map(r=>[r.gene]);
 sh.getRange('C80:C142').values=rows.map(()=>[cfg.unit]);
 sh.getRange('B80:B142').formulas=rows.map((r,i)=>[`=COUNT(D${sourceFirst+i}:${last}${sourceFirst+i})`]);
 sh.getRange(`D8:${last}70`).formulas=rows.map((r,i)=>species.map((sp,j)=>{
  const c=col(j+4),src=`${c}${sourceFirst+i}`,reference=`${c}$${sourceFirst}:${c}$${sourceLast}`;
  return `=IF(ISNUMBER(${src}),1+COUNTIFS(${reference},">"&${src})+(COUNTIFS(${reference},${src})-1)/2,"NA")`;
 }));
 sh.getRange('B8:B70').formulas=rows.map((r,i)=>[`=IF(COUNT(D${firstRow+i}:${last}${firstRow+i})=0,"NA",AVERAGE(D${firstRow+i}:${last}${firstRow+i}))`]);
 sh.getRange('C8:C70').formulas=rows.map((r,i)=>[`=COUNT(D${firstRow+i}:${last}${firstRow+i})`]);
 sh.getRange(`B8:${last}70`).format.horizontalAlignment='right';
 sh.getRange(`D8:${last}70`).setNumberFormat('0.0');sh.getRange('B8:B70').setNumberFormat('0.00');sh.getRange('C8:C70').setNumberFormat('0');
 sh.getRange(`D80:${last}142`).setNumberFormat(cfg.format);sh.getRange(`D80:${last}142`).format.horizontalAlignment='right';
 sh.getRange('B80:B142').setNumberFormat('0');
 for(let r=8;r<=70;r+=2)sh.getRange(`A${r}:C${r}`).format.fill=pale;
 for(let r=80;r<=142;r+=2)sh.getRange(`A${r}:${last}${r}`).format.fill='#F4F7FA';
 sh.getRange(`D8:${last}70`).conditionalFormats.add('colorScale',{
  colors:['#86B9D5','#DDECF4','#FFFFFF'],thresholds:[1,32,63]
 });
 for(const range of [`D8:${last}70`,`D80:${last}142`])sh.getRange(range).conditionalFormats.add('containsText',{
  text:'NA',format:{fill:'#ECEFF2',font:{color:'#79838C'}}
 });
 sh.getRange('A74').values=[['Underlying species-median values']];
 sh.getRange('A75').values=[[cfg.source]];
 sh.getRange('A76').values=[[cfg.definition]];
 sh.getRange('A77').values=[['The rank matrix above references these values. Missing values stay NA; genuine zeros remain zero.']];
 sh.getRange('A78').values=[['Editing values recalculates ranks and mean ranks. The saved row order is fixed; sort Mean rank ascending to reorder.']];
 sh.getRange(`A79:${last}79`).values=[['Gene / label','Species measured','Unit',...species.map(s=>s.common)]];
 for(const r of [7,79]){
  const h=sh.getRange(`A${r}:${last}${r}`);h.format={fill:headerColor,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:55};
  h.format.borders={insideVertical:{style:'thin',color:'#FFFFFF'}};
 }
 sh.getRange('A144').values=[['Calculation and matching notes']];sh.getRange('A144').format.font={name:'Arial',size:11,bold:true,color:headerColor};
 const notes=[
  'Case-insensitive exact gene symbols match the browser. No gene-alias substitution or verified orthology mapping.',
  'All 63 screenshot labels are retained, including noncoding RNA and pseudogene annotations.',
  cfg.key==='rna'?'One measurement per gene and sample was verified; species medians weight each sample equally.':'For multiple loci, use the highest species-median activity, then apply the region-availability rule.',
  cfg.key==='rna'?'The RNA sample BMR_Turk_Lu is excluded because it is lung.':'Selected loci with less than 50% of the requested 10 kb span inside the assembly are unavailable. Exactly 50% is retained.',
  cfg.key==='rna'?'RNA species Ellobius lutescens is matched to browser Bramus lutescens (Transcaucasian mole vole).':'No log transformation, gene-length adjustment, body-mass adjustment, or MLS filtering is used for these abundance ranks.',
  cfg.key==='rna'?'African spiny mouse: source species says Acomys cahirinus, but its common field says Acomys russatus; identity is unresolved.':'The locus TSS is the start on the plus strand or final body base on the minus strand, matching the activity panel.',
  'Mean rank is descriptive: each gene may be measured in a different set of species and ranked against a different number of labels.'
 ];
 sh.getRange(`A145:A${144+notes.length}`).values=notes.map(n=>[n]);
 sh.freezePanes.freezeRows(7);sh.freezePanes.freezeColumns(3);
}

wb.recalculate();
// Check all exported ranks and medians against the independent Python result.
for(const cfg of configs){
 const sh=wb.worksheets.getItem(cfg.name),rows=data[cfg.key];
 const ranks=sh.getRange(`D8:${last}70`).values,means=sh.getRange('B8:C70').values;
 for(let i=0;i<rows.length;i++){
  for(let j=0;j<species.length;j++){
   const want=rows[i].ranks[j]??'NA';if(ranks[i][j]!==want)throw Error(`${cfg.key} rank mismatch ${i},${j}: ${ranks[i][j]} vs ${want}`);
  }
  if(rows[i].mean_rank===null ? means[i][0]!=='NA' : Math.abs(means[i][0]-rows[i].mean_rank)>1e-10)throw Error('Mean rank mismatch');
  if(means[i][1]!==rows[i].n_species)throw Error('Coverage mismatch');
 }
 // Verify a change from NA to a very high source value updates its rank.
 const src=sh.getRange('D80'),old=src.values;
 src.values=[[1e12]];
 if(sh.getRange('D8').values[0][0]!==1)throw Error('Rank failed to update on a source-value change');
 src.values=old;
}
wb.recalculate();
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Formula error scan',maxChars:1000})).ndjson);
for(const cfg of configs){
 console.log((await wb.inspect({kind:'table',range:`'${cfg.name}'!A7:G11`,include:'values,formulas',tableMaxRows:5,tableMaxCols:7,maxChars:2200})).ndjson);
 for(const [name,range] of [['ranks','A1:J19'],['values','A74:J87']]){
  const preview=await wb.render({sheetName:cfg.name,range,scale:1.4,format:'png'});
  await fs.writeFile(path.join(out,`preview_${cfg.key}_${name}.png`),new Uint8Array(await preview.arrayBuffer()));
 }
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(out,'Gene_expression_and_TSS5kb_rank_matrices.xlsx'));
console.log('Saved workbook: two matrices; all 4,032 ranks and 126 mean ranks verified.');
