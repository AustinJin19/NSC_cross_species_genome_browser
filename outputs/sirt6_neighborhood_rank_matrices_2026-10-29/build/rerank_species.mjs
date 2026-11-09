import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {FileBlob, SpreadsheetFile} from '@oai/artifact-tool';

const source=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const out=path.resolve(source,'../sirt6_species_ranks_2026-10-30');
await fs.mkdir(out,{recursive:true});
const data=JSON.parse(await fs.readFile(path.join(source,'matrix_data.json'),'utf8'));
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(path.join(source,'Gene_expression_and_TSS5kb_rank_matrices.xlsx')));
const species=data.metadata.species;
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const configs=[{key:'rna',name:'RNA expression',title:'Species ranked by liver RNA expression',file:'RNA_expression_species_ranks'},
 {key:'atac',name:'TSS 5kb activity',title:'Species ranked by TSS ±5 kb gene activity',file:'TSS_5kb_activity_species_ranks'}];
const expected={};
for(const cfg of configs){
 const sh=wb.worksheets.getItem(cfg.name);
 const oldGenes=sh.getRange('A80:A142').values.map(r=>r[0]);
 const oldValues=sh.getRange('D80:AI142').values;
 const oldMap=new Map(oldGenes.map((g,i)=>[g,oldValues[i]]));
 const rows=[...data[cfg.key]].sort((a,b)=>a.gene.toUpperCase().localeCompare(b.gene.toUpperCase()));
 for(const row of rows){
  if(JSON.stringify(oldMap.get(row.gene))!==JSON.stringify(row.values.map(v=>v??'NA')))throw Error(`Source values differ: ${row.gene}`);
 }
 expected[cfg.key]=rows.map(row=>{
  const vals=row.values.filter(v=>v!==null);
  return {gene:row.gene,values:row.values,n_species:vals.length,
   highest_species:vals.length?species[row.values.indexOf(Math.max(...vals))].common:'NA',
   ranks:row.values.map(v=>v===null?'NA':1+vals.filter(x=>x>v).length+(vals.filter(x=>x===v).length-1)/2)};
 });
 sh.getRange('A2').values=[[cfg.title]];
 sh.getRange('A3').values=[['For each gene (row), rank species from highest to lowest. Rank 1 = highest. Exact ties use average ranks. NA = unavailable.']];
 sh.getRange('A5').values=[['Genes are alphabetical in both matrices. Highest species shows the first tied species when multiple species share the maximum.']];
 sh.getRange('B7').values=[['Highest species']];
 sh.getRange('B1:B157').format.columnWidth=29;
 sh.getRange('A8:A70').values=rows.map(r=>[r.gene]);
 sh.getRange('A80:A142').values=rows.map(r=>[r.gene]);
 sh.getRange('D80:AI142').values=rows.map(r=>r.values.map(v=>v??'NA'));
 sh.getRange('D8:AI70').formulas=rows.map((r,i)=>species.map((s,j)=>{
  const n=80+i,src=`${col(j+4)}${n}`,ref=`$D${n}:$AI${n}`;
  return `=IF(ISNUMBER(${src}),1+COUNTIFS(${ref},">"&${src})+(COUNTIFS(${ref},${src})-1)/2,"NA")`;
 }));
 sh.getRange('B8:B70').formulas=rows.map((r,i)=>{
  const n=80+i;return [`=IF(COUNT(D${n}:AI${n})=0,"NA",INDEX($D$79:$AI$79,1,MATCH(MAX(D${n}:AI${n}),D${n}:AI${n},0)))`];
 });
 sh.getRange('B8:B70').setNumberFormat('@');
 sh.getRange('B8:B70').format.horizontalAlignment='left';
 sh.getRange('A78').values=[['Editing the underlying values recalculates the species ranks across that gene’s row. Gene and species order stay fixed.']];
 sh.getRange('A151').values=[['Ranks compare species using the supplied signal normalization. CPM alone does not remove cross-species or batch effects.']];
 sh.getRange('A152').values=[['Ranks are recalculated across species from the same species-median values used in the 29 October 2026 workbook.']];
 const grid=sh.getRange('D8:AI70');grid.conditionalFormats.deleteAll();
 grid.conditionalFormats.add('colorScale',{colors:['#86B9D5','#DDECF4','#FFFFFF'],thresholds:[1,16.5,32]});
 grid.conditionalFormats.add('containsText',{text:'NA',format:{fill:'#ECEFF2',font:{color:'#79838C'}}});
 const tsv=[['Gene',...species.map(s=>s.common)],...expected[cfg.key].map(r=>[r.gene,...r.ranks])].map(r=>r.join('\t')).join('\n')+'\n';
 await fs.writeFile(path.join(out,cfg.file+'.tsv'),tsv);
}

wb.recalculate();
function verify(){
 for(const cfg of configs){
  const sh=wb.worksheets.getItem(cfg.name),rows=expected[cfg.key];
  const ranks=sh.getRange('D8:AI70').values,highest=sh.getRange('B8:C70').values,values=sh.getRange('D80:AI142').values;
  for(let i=0;i<rows.length;i++){
   if(JSON.stringify(ranks[i])!==JSON.stringify(rows[i].ranks))throw Error(`Rank mismatch: ${cfg.key} ${rows[i].gene}`);
   if(highest[i][0]!==rows[i].highest_species||highest[i][1]!==rows[i].n_species)throw Error(`Summary mismatch: ${cfg.key} ${rows[i].gene}`);
   if(JSON.stringify(values[i])!==JSON.stringify(rows[i].values.map(v=>v??'NA')))throw Error('Changed original median values');
  }
 }
}
verify();
for(const cfg of configs){
 const sh=wb.worksheets.getItem(cfg.name),row=sh.getRange('D80:AI80'),old=row.values;
 row.values=[[100,100,0,...Array(29).fill('NA')]];
 const test=sh.getRange('D8:G8').values[0];
 if(JSON.stringify(test)!==JSON.stringify([1.5,1.5,3,'NA']))throw Error(`Tie/zero/missing test failed: ${test}`);
 row.values=old;
}
wb.recalculate();verify();
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:10},summary:'Formula error scan',maxChars:1000})).ndjson);
for(const cfg of configs){
 console.log((await wb.inspect({kind:'table',range:`'${cfg.name}'!A7:H11`,include:'values,formulas',tableMaxRows:5,tableMaxCols:8,maxChars:1800})).ndjson);
 const preview=await wb.render({sheetName:cfg.name,range:'A1:J18',scale:1.4,format:'png'});
 await fs.writeFile(path.join(out,`preview_${cfg.key}.png`),new Uint8Array(await preview.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(out,'Species_ranks_RNA_and_TSS5kb.xlsx'));
await fs.writeFile(path.join(out,'species_rank_data.json'),JSON.stringify({species,rank_definition:'Rank species across each gene row, descending; average ties; NA excluded',source:source, ...expected},null,2)+'\n');
console.log('Saved two matrices: 63 genes × 32 species. Source values preserved. All ranks, ties, zeros and missing values verified.');
for(const cfg of configs){const r=expected[cfg.key].find(r=>r.gene==='SIRT6');console.log(cfg.key+' SIRT6 top species: '+r.highest_species);}
