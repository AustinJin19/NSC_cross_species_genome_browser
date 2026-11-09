import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const source=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const out=path.resolve(source,'../ordered_genes_by_species_2026-10-30');
await fs.mkdir(out,{recursive:true});
const data=JSON.parse(await fs.readFile(path.join(source,'matrix_data.json'),'utf8'));
const wb=Workbook.create(),species=data.metadata.species;
const col=n=>{let s='';for(;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const configs=[{key:'rna',name:'RNA expression',title:'Genes ranked within each species: liver RNA expression',unit:'Source log.tpm',format:'0.000'},
 {key:'atac',name:'TSS 5kb activity',title:'Genes ranked within each species: TSS ±5 kb activity',unit:'CPM·bp',format:'#,##0.0'}];
const verification=[],manifest={species,source:'Species median values from 29 October 2026; no change to measurements.',ranking:'Highest to lowest within each species. Exact ties ordered alphabetically. Missing measurements are not ranked.'};
for(const cfg of configs){
 const sh=wb.worksheets.add(cfg.name),genes=[...data[cfg.key]].sort((a,b)=>a.gene.toUpperCase().localeCompare(b.gene.toUpperCase()));
 const lists=species.map((sp,j)=>genes.filter(g=>g.values[j]!==null).sort((a,b)=>b.values[j]-a.values[j]||a.gene.toUpperCase().localeCompare(b.gene.toUpperCase())));
 const length=Math.max(...lists.map(l=>l.length)),listFirst=9,listLast=8+length;
 const sourceTitle=listLast+3,sourceHeader=sourceTitle+3,sourceFirst=sourceHeader+1,sourceLast=sourceFirst+62;
 const rankTitle=sourceLast+3,rankHeader=rankTitle+2,rankFirst=rankHeader+1,rankLast=rankFirst+62;
 sh.showGridLines=false;sh.tabColor=cfg.key==='rna'?'#183F62':'#375F49';
 const all=sh.getRange(`A1:AG${rankLast+3}`);all.format.font={name:'Arial',size:10,color:'#203246'};all.format.rowHeight=20;all.format.verticalAlignment='center';
 sh.getRange(`A1:A${rankLast+3}`).format.columnWidth=22;
 sh.getRange(`B1:AG${rankLast+3}`).format.columnWidth=19;
 sh.getRange('A2').values=[[cfg.title]];sh.getRange('A2:AG2').format.font={name:'Arial',size:14,bold:true,color:'#183F62'};
 sh.getRange('A3').values=[['Each species column lists gene names from highest (1) to lowest, ranked among the 63 screenshot labels with measurements.']];
 sh.getRange('A4').values=[['Human and bowhead whale are excluded. Macaque is retained; its RNA list is empty because matching liver RNA is unavailable.']];
 sh.getRange('A5').values=[['Equal values are ordered alphabetically. Empty list cells mean no additional measured genes. Source NA values are not ranked.']];
 sh.getRange('A6').values=[[cfg.key==='rna'?'Ranked by the species median of source liver log.tpm values.':'Ranked by species-median integrated CPM in TSS ±5 kb (10 kb total).']];
 sh.getRange('A7:AG7').values=[['Rank',...species.map(s=>s.common)]];
 sh.getRange('A8').values=[['Genes ranked']];
 sh.getRange('B8:AG8').formulas=[species.map((s,j)=>`=COUNT(${col(j+2)}${sourceFirst}:${col(j+2)}${sourceLast})`)];
 sh.getRange('A8:AG8').format.fill='#EDF3F8';sh.getRange('A8:AG8').format.font={name:'Arial',size:10,bold:true,color:'#536579'};
 sh.getRange(`A${listFirst}:A${listLast}`).values=Array.from({length},(_,i)=>[i+1]);
 sh.getRange(`B${listFirst}:AG${listLast}`).formulas=Array.from({length},(_,i)=>species.map((s,j)=>{
  const c=col(j+2),row=listFirst+i;
  return `=IF($A${row}<=${c}$8,INDEX($A$${sourceFirst}:$A$${sourceLast},MATCH($A${row},${c}$${rankFirst}:${c}$${rankLast},0)),"")`;
 }));
 sh.getRange(`A${sourceTitle}`).values=[['Underlying species-median values']];
 sh.getRange(`A${sourceTitle+1}`).values=[[cfg.key==='rna'?'Sources: CrossSpecies_Count.xlsx (liver), rnaLongTPM.xlsx (mouse liver supplement).':'Sources: browser gene-locus index and sample CPM bigWigs. Quantified 29 October 2026.']];
 sh.getRange(`A${sourceTitle+2}`).values=[[cfg.key==='rna'?'Units: source log.tpm; no additional logarithm. Original gene-symbol matching and source exclusions are retained.':'Units: CPM·bp. Highest-median locus selected first; available window must be at least 50%. No body-mass or length correction.']];
 sh.getRange(`A${sourceHeader}:AG${sourceHeader}`).values=[['Gene / label',...species.map(s=>s.common)]];
 sh.getRange(`A${sourceFirst}:A${sourceLast}`).values=genes.map(r=>[r.gene]);
 sh.getRange(`B${sourceFirst}:AG${sourceLast}`).values=genes.map(r=>r.values.map(v=>v??'NA'));
 sh.getRange(`B${sourceFirst}:AG${sourceLast}`).setNumberFormat(cfg.format);
 sh.getRange(`A${rankTitle}`).values=[['Order positions used for each species list']];
 sh.getRange(`A${rankTitle+1}`).values=[['Descending values. Alphabetical gene order breaks exact ties. Editing source values updates both order positions and gene lists.']];
 sh.getRange(`A${rankHeader}:AG${rankHeader}`).values=[['Gene / label',...species.map(s=>s.common)]];
 sh.getRange(`A${rankFirst}:A${rankLast}`).values=genes.map(r=>[r.gene]);
 sh.getRange(`B${rankFirst}:AG${rankLast}`).formulas=genes.map((r,i)=>species.map((s,j)=>{
  const c=col(j+2),cell=`${c}${sourceFirst+i}`;
  return `=IF(ISNUMBER(${cell}),COUNTIFS(${c}$${sourceFirst}:${c}$${sourceLast},">"&${cell})+COUNTIFS(${c}$${sourceFirst}:${cell},${cell}),"NA")`;
 }));
 for(const row of [7,sourceHeader,rankHeader])sh.getRange(`A${row}:AG${row}`).format={fill:'#183F62',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:55};
 for(const row of [sourceTitle,rankTitle])sh.getRange(`A${row}:AG${row}`).format.font={name:'Arial',size:12,bold:true,color:'#183F62'};
 for(const [a,b] of [[listFirst,listLast],[sourceFirst,sourceLast],[rankFirst,rankLast]])for(let row=a;row<=b;row+=2)sh.getRange(`A${row}:AG${row}`).format.fill='#F4F7FA';
 sh.getRange(`B${sourceFirst}:AG${sourceLast}`).conditionalFormats.add('containsText',{text:'NA',format:{fill:'#ECEFF2',font:{color:'#79838C'}}});
 sh.getRange(`B${rankFirst}:AG${rankLast}`).conditionalFormats.add('containsText',{text:'NA',format:{fill:'#ECEFF2',font:{color:'#79838C'}}});
 sh.getRange(`B${listFirst}:AG${listLast}`).conditionalFormats.add('cellIs',{operator:'equal',formula:'"SIRT6"',format:{fill:'#D9EAF5',font:{bold:true,color:'#183F62'}}});
 sh.getRange(`A${rankLast+2}`).values=[['Source caveat: African spiny mouse RNA species field says Acomys cahirinus, but its common-name field says Acomys russatus.']];
 sh.getRange(`A${rankLast+3}`).values=[['Matching uses the original symbols from the previous workbook. The new browser ortholog mapping has not been applied to these measurements.']];
 sh.freezePanes.freezeRows(8);sh.freezePanes.freezeColumns(1);
 const expected=Array.from({length},(_,i)=>lists.map(l=>l[i]?.gene??''));
 verification.push({sh,expected,sourceFirst,sourceLast,listFirst,listLast,genes});
 const rows=[['Rank',...species.map(s=>s.common)],...expected.map((r,i)=>[i+1,...r])];
 await fs.writeFile(path.join(out,cfg.key+'_ordered_gene_lists.tsv'),rows.map(r=>r.join('\t')).join('\n')+'\n');
 manifest[cfg.key]=species.map((sp,j)=>({species:sp.common,genes:lists[j].map(g=>g.gene)}));
}
wb.recalculate();
function check(){for(const v of verification){if(JSON.stringify(v.sh.getRange(`B${v.listFirst}:AG${v.listLast}`).values)!==JSON.stringify(v.expected))throw Error('Ordered gene list does not match sorted original measurements');if(JSON.stringify(v.sh.getRange(`B${v.sourceFirst}:AG${v.sourceLast}`).values)!==JSON.stringify(v.genes.map(r=>r.values.map(x=>x??'NA'))))throw Error('Original values changed');}}
check();
for(const v of verification){const cell=v.sh.getRange(`B${v.sourceFirst}`),old=cell.values;cell.values=[[1e12]];if(v.sh.getRange('B9').values[0][0]!==v.genes[0].gene)throw Error('List did not reorder when value changed');cell.values=old;}
wb.recalculate();check();
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:10},maxChars:800})).ndjson);
for(const cfg of configs){const preview=await wb.render({sheetName:cfg.name,range:'A1:H20',scale:1.4,format:'png'});await fs.writeFile(path.join(out,`preview_${cfg.key}.png`),new Uint8Array(await preview.arrayBuffer()));}
await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(out,'Ordered_gene_lists_by_species.xlsx'));
await fs.writeFile(path.join(out,'ordered_lists.json'),JSON.stringify(manifest,null,2)+'\n');
console.log('Saved two ordered gene-name matrices; both checked against original values and tested for live reordering.');
